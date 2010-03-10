/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>

/* fftw header */
#include <fftw3.h>

/* own header */
#include "poisson.h"


/*>>>>>>>>>>>>> Make matrices for solving the particular part of the radial equation <<<<<<<<<<<<<<<*/

/*************************************************/
/* radial[z][L] is a pointer to an nr*nr matrix  */
/*************************************************/
gsl_matrix ***radial_matrix_alloc(int nz, int nr, int nt)
{
  int z;
  int L;
  gsl_matrix ***radial;
  
  /* allocate memory for the array of nz pointers that each point to an array of nt pointers */
  radial = (gsl_matrix ***) malloc(nz*sizeof(gsl_matrix **));
  
  for (z = 0; z < nz; z++) {
    /* allocate memory for the array of nt pointers */
    radial[z] = (gsl_matrix **) malloc(nt*sizeof(gsl_matrix *));
    for (L = 0; L < nt; L++) {
      /* allocate memory for each nr*nr matrix */
      radial[z][L] = gsl_matrix_alloc(nr, nr);
    }
  }
  
  return radial;
}



void radial_matrix_set(int nz, int nt, gsl_vector *alpha_v, gsl_vector *beta_v, gsl_matrix ***radial)
{
  int z;
  int L;
 
  /* set matrices in kernel for even L */
  for (L = 0; L < nt; L += 2) {
    set_A_kernel_even(radial[0][L], L);
  }
  /* set matrices in kernel for odd L */
  for (L = 1; L < nt; L += 2) {
    set_A_kernel_odd(radial[0][L], L);
  }
  /* set matrices in shells */
  for (z = 1; z <= nz-2; z++) { /* only evaluated if nz >= 3 */
    for (L = 0; L < nt; L++) {
      set_A_shell(radial[z][L], L, gsl_vector_get(alpha_v, z), gsl_vector_get(beta_v, z));
    }
  }
  /* set matrices in external domain */
  for (L = 0; L < nt; L++) {
    set_A_ext(radial[nz-1][L], L);
  }
}

void radial_matrix_free(int nz, int nt, gsl_matrix ***radial)
{
  int z;
  int L;

  for (z = 0; z < nz; z++) {
    for (L = 0; L < nt; L++) {
      /* free memory for each nr*nr matrix */
      gsl_matrix_free(radial[z][L]);
    }
    /* free memory for the array of nt pointers */
    free(radial[z]);
  }
  free(radial);
}



/*******************************************************************/
/* Effective source (including residual) for Poisson equation.     */
/* Eq. 111 of Bonazzola, Gourgoulhon, Marck PRD 58, 104020 (1998). */
/*******************************************************************/
void effective_source(scalar3d *s_eff_scalar3d, scalar3d *f_scalar3d, scalar3d *s_scalar3d, scalar3d *s_eff_jm1_scalar3d, scalar3d *s_eff_jm2_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k;
  int nz, nr, nt, np;
  
  scalar3d *j1_scalar3d;
  scalar3d *j2_scalar3d;
  scalar3d *j3_scalar3d;
  scalar3d *residual_scalar3d;
  scalar3d *a_scalar3d;
  gsl_vector *maxa_l_vector;

  double j1, j2, j3;
  double alpha, a, maxa_l;
  double s, residual_d, s_eff_jm1, s_eff_jm2, s_eff;
  
  double lambda = 0.5; /* relaxation parameter */

  nz = s_eff_scalar3d->nz;
  nr = s_eff_scalar3d->nr;
  nt = s_eff_scalar3d->nt;
  np = s_eff_scalar3d->np;
  
  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j3_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  residual_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  a_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  maxa_l_vector = gsl_vector_alloc(nz);

  jacobian1(j1_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  jacobian2(j2_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  jacobian3(j3_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  residual(residual_scalar3d, f_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  /*print_scalar3d(residual_scalar3d);*/

  /* evaluate a and maxa_l */
  for ( z = 0; z < nz; z++ ) {
    alpha = gsl_vector_get(alpha_vector, z);
    maxa_l = 0.0;
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	  j2 = scalar3d_get(j2_scalar3d, z, i, j, k);
	  j3 = scalar3d_get(j3_scalar3d, z, i, j, k);
	  
	  a = alpha*alpha * (1.0 + j2*j2 + j3*j3) / (j1*j1);
	  scalar3d_set(a_scalar3d, z, i, j, k, a);
	  maxa_l = MAX(a, maxa_l);
	}
      }
    }
    gsl_vector_set(maxa_l_vector, z, maxa_l);
  }
/*   print_scalar3d(a_scalar3d); */
/*   print_vector(maxa_l_vector); */

  /* set new effective source */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  s = scalar3d_get(s_scalar3d, z, i, j, k);
	  residual_d = scalar3d_get(residual_scalar3d, z, i, j, k);
	  a = scalar3d_get(a_scalar3d, z, i, j, k);
	  s_eff_jm1 = scalar3d_get(s_eff_jm1_scalar3d, z, i, j, k);
	  s_eff_jm2 = scalar3d_get(s_eff_jm2_scalar3d, z, i, j, k);
	  maxa_l = gsl_vector_get(maxa_l_vector, z);

	  s_eff = ( s + residual_d + (maxa_l - a)*(lambda*s_eff_jm1 + (1.0 - lambda)*s_eff_jm2) ) / maxa_l;
	  scalar3d_set(s_eff_scalar3d, z, i, j, k, s_eff);
	}
      }
    }
  }

  scalar3d_free(j1_scalar3d);
  scalar3d_free(j2_scalar3d);
  scalar3d_free(j3_scalar3d);
  scalar3d_free(residual_scalar3d);
  scalar3d_free(a_scalar3d);
  gsl_vector_free(maxa_l_vector);
}


/* /\*******************************************************************\/ */
/* /\* Evaluate homogeneous part of Poisson Eq. at an arbitrary point. *\/ */
/* /\*******************************************************************\/ */
/* void evaluate_homogeneous(scalar3d *homo_scalar3d, ylm_coeff *homo_grow_ylm_coeff, ylm_coeff *homo_decay_ylm_coeff, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d) */
/* { */
/*   int z, i, j, k; */
/*   int L, m, imag; */
/*   int nz, nr, nt, np; */

/*   scalar3d *r_scalar3d; */
/*   double roru_ijk; */
/*   double theta_j, x_j, phi_k; */
/*   double fhomo; */
/*   double phipart; */
/*   double ylm_jk; */
/*   double A_lm, B_lm; */

/*   nz = homo_scalar3d->nz; */
/*   nr = homo_scalar3d->nr; */
/*   nt = homo_scalar3d->nt; */
/*   np = homo_scalar3d->np; */

/*   r_scalar3d = scalar3d_alloc(nz, nr, nt, np); */

/*   /\* determine radius or 1/radius for each gridpoint *\/ */
/*   rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */
  
/* /\*   print_ylm_coeff(homo_grow_ylm_coeff); *\/ */
/* /\*   print_ylm_coeff(homo_decay_ylm_coeff); *\/ */

/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  roru_ijk = scalar3d_get(r_scalar3d, z, i, j, k); */
/* 	  theta_j = PI*j/(nt-1);  */
/* 	  x_j = cos(theta_j); */
/* 	  phi_k = 2*PI*k/np;	 */
	  
/* 	  /\* begin sum over all spherical harmonics *\/	   */
/* 	  fhomo = 0.0; */
/* 	  for ( imag=0; imag<=1; imag++ ) { */
/* 	    for ( m = imag; m < nt-imag; m++ ) { */
/* 	      if(imag==0) { */
/* 		phipart = cos(m*phi_k); */
/* 	      } else { */
/* 		phipart = sin(m*phi_k); */
/* 	      } */
/* 	      for ( L = m; L < nt-m%2; L++ ) { */
/* 		A_lm = ylm_coeff_get(homo_grow_ylm_coeff, z, 0, L, m, imag); */
/* 		B_lm = ylm_coeff_get(homo_decay_ylm_coeff, z, 0, L, m, imag); */
/* 		ylm_jk = gsl_sf_legendre_sphPlm(L, m, x_j) * phipart; */
/* 		if(z==0) { /\* in kernel *\/ */
/* 		  fhomo += A_lm * pow(roru_ijk, L) * ylm_jk; */
/* 		} else if(z==nz-1) { /\* in external zone *\/ */
/* 		  fhomo += B_lm * pow(roru_ijk, (L+1)) * ylm_jk;		  */
/* 		} else { /\* in a shell *\/ */
/* 		  fhomo += (A_lm * pow(roru_ijk, L) + B_lm * pow(roru_ijk, -(L+1))) * ylm_jk;		  */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	  scalar3d_set(homo_scalar3d, z, i, j, k, fhomo); */

/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   r_scalar3d = scalar3d_alloc(nz, nr, nt, np); */
/* } */


/* /\*************************************************************************************************\/ */
/* /\* Solve one iteration of the Poisson equation for a source evaluated at grid points. *\/ */
/* /\*************************************************************************************************\/ */
/* void poisson_iteration(scalar3d *field_scalar3d, scalar3d *s_eff_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d) */
/* { */
/*   int z, i, j, k; */
/*   int L, m, imag; */
/*   int nz, nr, nt, np; */

/*   gsl_matrix ***radial_matrix; */
/*   gsl_matrix **fouriertoylm_matrix; */
/*   gsl_matrix **ylmtofourier_matrix; */

/*   coeff *source_coeff; */
/*   ylm_coeff *source_ylm_coeff; */

/*   gsl_matrix *source_matrix; */
/*   gsl_matrix *particular_matrix; */
/*   ylm_coeff *particular_ylm_coeff; */
/*   coeff *particular_coeff; */
/*   scalar3d *particular_scalar3d; */

/*   gsl_vector *homo_grow_vector; */
/*   gsl_vector *homo_decay_vector; */
/*   ylm_coeff *homo_grow_ylm_coeff; */
/*   ylm_coeff *homo_decay_ylm_coeff; */
/*   scalar3d *homo_scalar3d; */

/*   nz = field_scalar3d->nz; */
/*   nr = field_scalar3d->nr; */
/*   nt = field_scalar3d->nt; */
/*   np = field_scalar3d->np; */
  
/*   /\* allocate and make matrices for solving the radial poisson equation *\/ */
/*   radial_matrix = radial_matrix_alloc(nz, nr, nt); */
/*   radial_matrix_set(nz, nt, alpha_vector, beta_vector, radial_matrix); */
  
/*   /\* allocate and make matrices for fourier <--> spherical harmonic transforms *\/ */
/*   fouriertoylm_matrix = fouriertoylm_matrix_alloc(nt); */
/*   ylmtofourier_matrix = ylmtofourier_matrix_alloc(nt); */
/*   fouriertoylm_matrix_set(fouriertoylm_matrix); */
/*   ylmtofourier_matrix_set(ylmtofourier_matrix); */
  
/*   /\* allocate other things *\/ */
/*   source_coeff = coeff_alloc(nz, nr, nt, np); */
/*   source_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np); */
/*   source_matrix = gsl_matrix_alloc(nz, nr); */
/*   particular_matrix = gsl_matrix_alloc(nz, nr); */
/*   particular_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np); */
/*   particular_coeff = coeff_alloc(nz, nr, nt, np); */
/*   particular_scalar3d = scalar3d_alloc(nz, nr, nt, np); */
/*   homo_grow_vector = gsl_vector_alloc(nz); */
/*   homo_decay_vector = gsl_vector_alloc(nz); */
/*   homo_grow_ylm_coeff = ylm_coeff_alloc(nz, 1, nt, np); */
/*   homo_decay_ylm_coeff = ylm_coeff_alloc(nz, 1, nt, np); */
/*   homo_scalar3d = scalar3d_alloc(nz, nr, nt, np); */

/*   /\* go to fourier series coefficients *\/ */
/*   gridtofourier(source_coeff, s_eff_scalar3d, 0, 0); */
    
/*   /\* go from fourier series to spherical harmonics *\/ */
/*   transform_fouriertoylm(source_coeff, source_ylm_coeff, fouriertoylm_matrix); */

/*   /\*>>>>>>>>>>>>>>>> Solve particular and homogeneous parts of poisson equation. <<<<<<<<<<<<<<<<<<<<*\/ */
/*   /\* Particular solution is returned as coefficients of a_iLm*T_i(xi)*Y_lm(theta, phi) basis.        *\/ */
/*   /\* Homogeneous solution is returned as coefficients of (A_Lm*r^L + B_Lm*r^(-L-1))*Y_lm(theta, phi) *\/ */
/*   for ( imag=0; imag<=1; imag++ ) { */
/*     for ( m = imag; m < nt-imag; m++ ) { */
/*       for ( L = m; L < nt-m%2; L++ ) { */
/* 	/\* Solve particular part *\/ */
/* 	for ( z = 0; z < nz; z++ ) { */
/* 	  for ( i = 0; i < nr; i++ ) { */
/* 	    gsl_matrix_set(source_matrix, z, i, ylm_coeff_get(source_ylm_coeff, z, i, L, m, imag)); */
/* 	  } */
/* 	} */
/* 	gsl_matrix_set_zero(particular_matrix); */
/* 	solve_radial_particular(L, radial_matrix, alpha_vector, beta_vector, source_matrix, particular_matrix); */
/* 	for ( z = 0; z < nz; z++ ) { */
/* 	  for ( i = 0; i < nr; i++ ) { */
/* 	    ylm_coeff_set(particular_ylm_coeff, z, i, L, m, imag, gsl_matrix_get(particular_matrix, z, i)); */
/* 	  } */
/* 	} */
/* 	/\* Solve homogeneous part *\/ */
/* 	solve_radial_homogeneous(L, alpha_vector, beta_vector, particular_matrix, homo_grow_vector, homo_decay_vector); */
/* 	for ( z = 0; z < nz; z++ ) { */
/* 	  ylm_coeff_set(homo_grow_ylm_coeff, z, 0, L, m, imag, gsl_vector_get(homo_grow_vector, z)); */
/* 	  ylm_coeff_set(homo_decay_ylm_coeff, z, 0, L, m, imag, gsl_vector_get(homo_decay_vector, z)); */
/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   /\* evaluate particular solution on grid *\/ */
/*   transform_ylmtofourier(particular_ylm_coeff, particular_coeff, ylmtofourier_matrix); */
/*   fouriertogrid(particular_scalar3d, particular_coeff, 0, 0); */

/*   /\* evaluate homogeneous solution on grid *\/ */
/*   evaluate_homogeneous(homo_scalar3d, homo_grow_ylm_coeff, homo_decay_ylm_coeff, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */

/*   /\* field = particular + homogeneous *\/ */
/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  scalar3d_set(field_scalar3d, z, i, j, k,  */
/* 		       scalar3d_get(particular_scalar3d, z, i, j, k) + scalar3d_get(homo_scalar3d, z, i, j, k)); */
/* 	} */
/*       } */
/*     } */
/*   }   */
  
/*   /\* free memory *\/ */
/*   radial_matrix_free(nz, nt, radial_matrix); */
/*   fouriertoylm_matrix_free(fouriertoylm_matrix); */
/*   ylmtofourier_matrix_free(ylmtofourier_matrix); */
/*   coeff_free(source_coeff); */
/*   ylm_coeff_free(source_ylm_coeff); */
/*   gsl_matrix_free(source_matrix); */
/*   gsl_matrix_free(particular_matrix); */
/*   ylm_coeff_free(particular_ylm_coeff); */
/*   coeff_free(particular_coeff); */
/*   scalar3d_free(particular_scalar3d); */
/*   gsl_vector_free(homo_grow_vector); */
/*   gsl_vector_free(homo_decay_vector); */
/*   ylm_coeff_free(homo_grow_ylm_coeff); */
/*   ylm_coeff_free(homo_decay_ylm_coeff); */
/*   scalar3d_free(homo_scalar3d); */
/* } */


/*************************************************************************************************/
/* Solve one iteration of the Poisson equation for a source evaluated at grid points. */
/*************************************************************************************************/
void poisson_iteration(scalar3d *field_scalar3d, scalar3d *s_eff_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k;
  int L, m, imag;
  int nz, nr, nt, np;

  gsl_matrix ***radial_matrix;
  gsl_matrix **fouriertoylm_matrix;
  gsl_matrix **ylmtofourier_matrix;

  coeff *source_coeff;
  ylm_coeff *source_ylm_coeff;

  gsl_matrix *source_matrix;
  gsl_matrix *particular_matrix;
  ylm_coeff *particular_ylm_coeff;

  gsl_vector *homo_grow_vector;
  gsl_vector *homo_decay_vector;
  ylm_coeff *homo_grow_ylm_coeff;
  ylm_coeff *homo_decay_ylm_coeff;
  ylm_coeff *homogeneous_ylm_coeff;

  ylm_coeff *field_ylm_coeff;
  coeff *field_coeff;

  nz = field_scalar3d->nz;
  nr = field_scalar3d->nr;
  nt = field_scalar3d->nt;
  np = field_scalar3d->np;
  
  /* allocate and make matrices for solving the radial poisson equation */
  radial_matrix = radial_matrix_alloc(nz, nr, nt);
  radial_matrix_set(nz, nt, alpha_vector, beta_vector, radial_matrix);
  
  /* allocate and make matrices for fourier <--> spherical harmonic transforms */
  fouriertoylm_matrix = fouriertoylm_matrix_alloc(nt);
  ylmtofourier_matrix = ylmtofourier_matrix_alloc(nt);
  fouriertoylm_matrix_set(fouriertoylm_matrix);
  ylmtofourier_matrix_set(ylmtofourier_matrix);
  
  /* allocate other things */
  source_coeff = coeff_alloc(nz, nr, nt, np);
  source_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  source_matrix = gsl_matrix_alloc(nz, nr);
  particular_matrix = gsl_matrix_alloc(nz, nr);
  particular_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  homo_grow_vector = gsl_vector_alloc(nz);
  homo_decay_vector = gsl_vector_alloc(nz);
  homo_grow_ylm_coeff = ylm_coeff_alloc(nz, 1, nt, np);
  homo_decay_ylm_coeff = ylm_coeff_alloc(nz, 1, nt, np);
  homogeneous_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  field_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);


  /* go to fourier series coefficients */
  gridtofourier(source_coeff, s_eff_scalar3d, 0, 0);
 /*  print_coeff(source_coeff);  */
    
  /* go from fourier series to spherical harmonics */
  transform_fouriertoylm(source_coeff, source_ylm_coeff, fouriertoylm_matrix);

/*   print_ylm_coeff(source_ylm_coeff); */

  /*>>>>>>>>>>>>>>>> Solve particular and homogeneous parts of poisson equation. <<<<<<<<<<<<<<<<<<<<*/
  /* Particular solution is returned as coefficients of a_iLm*T_i(xi)*Y_lm(theta, phi) basis.        */
  /* Homogeneous solution is returned as coefficients of (A_Lm*r^L + B_Lm*r^(-L-1))*Y_lm(theta, phi) */
  for ( imag=0; imag<=1; imag++ ) {
    for ( m = imag; m < nt-imag; m++ ) {
      for ( L = m; L < nt-m%2; L++ ) {
	/* Solve particular part */
	for ( z = 0; z < nz; z++ ) {
	  for ( i = 0; i < nr; i++ ) {
	    gsl_matrix_set(source_matrix, z, i, ylm_coeff_get(source_ylm_coeff, z, i, L, m, imag));
	  }
	}
	gsl_matrix_set_zero(particular_matrix);
	solve_radial_particular(L, radial_matrix, alpha_vector, beta_vector, source_matrix, particular_matrix);
	for ( z = 0; z < nz; z++ ) {
	  for ( i = 0; i < nr; i++ ) {
	    ylm_coeff_set(particular_ylm_coeff, z, i, L, m, imag, gsl_matrix_get(particular_matrix, z, i));
	  }
	}
	/* Solve homogeneous part */
	solve_radial_homogeneous(L, alpha_vector, beta_vector, particular_matrix, homo_grow_vector, homo_decay_vector);
	for ( z = 0; z < nz; z++ ) {
	  ylm_coeff_set(homo_grow_ylm_coeff, z, 0, L, m, imag, gsl_vector_get(homo_grow_vector, z));
	  ylm_coeff_set(homo_decay_ylm_coeff, z, 0, L, m, imag, gsl_vector_get(homo_decay_vector, z));
	}
      }
    }
  }

  /* Convert homogeneous solution from radial functions to coefficients in (xi, theta, phi) coordinate system */  
  homogeneoustochebyshev(homo_grow_ylm_coeff, homo_decay_ylm_coeff, alpha_vector, beta_vector, homogeneous_ylm_coeff);
  
  /* field = particular + homogeneous */
  for ( imag=0; imag<=1; imag++ ) {
    for ( m = imag; m < nt-imag; m++ ) {
      for ( L = m; L < nt-m%2; L++ ) {
	for ( z = 0; z < nz; z++ ) {
	  for ( i = 0; i < nr; i++ ) {
	    ylm_coeff_set(field_ylm_coeff, z, i, L, m, imag,
			  ylm_coeff_get(particular_ylm_coeff, z, i, L, m, imag) + ylm_coeff_get(homogeneous_ylm_coeff, z, i, L, m, imag));
	  }
	}
      }
    }
  }


  /* evaluate solution on grid */
  transform_ylmtofourier(field_ylm_coeff, field_coeff, ylmtofourier_matrix);
  fouriertogrid(field_scalar3d, field_coeff, 0, 0);

  
  /* free memory */
  radial_matrix_free(nz, nt, radial_matrix);
  fouriertoylm_matrix_free(fouriertoylm_matrix);
  ylmtofourier_matrix_free(ylmtofourier_matrix);
  coeff_free(source_coeff);
  ylm_coeff_free(source_ylm_coeff);
  gsl_matrix_free(source_matrix);
  gsl_matrix_free(particular_matrix);
  ylm_coeff_free(particular_ylm_coeff);
  gsl_vector_free(homo_grow_vector);
  gsl_vector_free(homo_decay_vector);
  ylm_coeff_free(homo_grow_ylm_coeff);
  ylm_coeff_free(homo_decay_ylm_coeff);
  ylm_coeff_free(homogeneous_ylm_coeff);
  ylm_coeff_free(field_ylm_coeff);
  coeff_free(field_coeff);
}


/* /\*************************************************************************************************\/ */
/* /\* Solve the Poisson equation for a source already decomposed into T_i(xi) and Y_l^m(theta, phi) *\/ */
/* /\*************************************************************************************************\/ */
/* void solve_poisson_spherical(ylm_coeff *field_ylm_coeff, ylm_coeff *source_ylm_coeff, gsl_vector *alpha_v, gsl_vector *beta_v) */
/* { */
/*   int imag, z, i, L, m; */
/*   int nz, nr, nt, np; */
/*   gsl_vector *homo_grow_v; */
/*   gsl_vector *homo_decay_v; */
/*   gsl_matrix *source_matrix; */
/*   gsl_matrix *particular_matrix; */
/*   gsl_matrix ***radial_matrix; */
/*   ylm_coeff *particular_ylm_coeff; */
/*   ylm_coeff *homogeneous_ylm_coeff; */
/*   ylm_coeff *homo_grow_ylm_coeff; */
/*   ylm_coeff *homo_decay_ylm_coeff; */

/*   nz = field_ylm_coeff->nz; */
/*   nr = field_ylm_coeff->nr; */
/*   nt = field_ylm_coeff->nt; */
/*   np = field_ylm_coeff->np; */

/*   homo_grow_v = gsl_vector_alloc(nz); */
/*   homo_decay_v = gsl_vector_alloc(nz); */
/*   source_matrix = gsl_matrix_alloc(nz, nr); */
/*   particular_matrix = gsl_matrix_alloc(nz, nr); */

/*   particular_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np); */
/*   homogeneous_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np); */
/*   homo_grow_ylm_coeff = ylm_coeff_alloc(nz, 1, nt, np); */
/*   homo_decay_ylm_coeff = ylm_coeff_alloc(nz, 1, nt, np); */
  
/*   /\* make matrices for solving the radial poisson equation *\/ */
/*   radial_matrix = radial_matrix_alloc(nz, nr, nt); */
/*   radial_matrix_set(nz, nt, alpha_v, beta_v, radial_matrix); */
  
/*   for ( imag=0; imag<=1; imag++ ) { */
/*     for ( m = imag; m < nt-imag; m++ ) { */
/*       for ( L = m; L < nt-m%2; L++ ) { */
/* 	/\* Solve each radial equation *\/ */
/* 	for ( z = 0; z < nz; z++ ) { */
/* 	  for ( i = 0; i < nr; i++ ) { */
/* 	    gsl_matrix_set(source_matrix, z, i, ylm_coeff_get(source_ylm_coeff, z, i, L, m, imag)); */
/* 	  } */
/* 	} */
/* 	/\*printf("imag=%d\tm=%d\tL=%d:\n", imag, m, L); */
/* 	  print_matrix(source_matrix);*\/ */
/* 	/\* Solve particular part *\/ */
/* 	/\*printf("imag=%d\tm=%d\tL=%d:\n", imag, m, L);*\/ */
/* 	gsl_matrix_set_zero(particular_matrix); */
/* 	solve_radial_particular(L, radial_matrix, alpha_v, beta_v, source_matrix, particular_matrix); */
/* 	for ( z = 0; z < nz; z++ ) { */
/* 	  for ( i = 0; i < nr; i++ ) { */
/* 	    ylm_coeff_set(particular_ylm_coeff, z, i, L, m, imag, gsl_matrix_get(particular_matrix, z, i)); */
/* 	  } */
/* 	} */
/* 	/\* Solve homogeneous part *\/ */
/* 	solve_radial_homogeneous(L, alpha_v, beta_v, particular_matrix, homo_grow_v, homo_decay_v); */
/* 	for ( z = 0; z < nz; z++ ) { */
/* 	  ylm_coeff_set(homo_grow_ylm_coeff, z, 0, L, m, imag, gsl_vector_get(homo_grow_v, z)); */
/* 	  ylm_coeff_set(homo_decay_ylm_coeff, z, 0, L, m, imag, gsl_vector_get(homo_decay_v, z)); */
/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   /\*print_ylm_coeff(homo_grow_ylm_coeff); */
/*     print_ylm_coeff(homo_decay_ylm_coeff);*\/ */
/*   /\* THIS CAUSES A BUS ERROR!!!!!!!!!!!!!! *\/ */
/*   /\* Convert homogeneous solution from radial functions to coefficients in (xi, theta, phi) coordinate system *\/ */
/*   homogeneoustochebyshev(homo_grow_ylm_coeff, homo_decay_ylm_coeff, alpha_v, beta_v, homogeneous_ylm_coeff); */
/*   /\*printf("here\n");*\/ */
/*   /\*print_ylm_coeff(homogeneous_ylm_coeff);*\/ */

/*   /\* field = particular + homogeneous *\/ */
/*   for ( imag=0; imag<=1; imag++ ) { */
/*     for ( m = imag; m < nt-imag; m++ ) { */
/*       for ( L = m; L < nt-m%2; L++ ) { */
/* 	for ( z = 0; z < nz; z++ ) { */
/* 	  for ( i = 0; i < nr; i++ ) { */
/* 	    ylm_coeff_set(field_ylm_coeff, z, i, L, m, imag, */
/* 			  ylm_coeff_get(particular_ylm_coeff, z, i, L, m, imag) + ylm_coeff_get(homogeneous_ylm_coeff, z, i, L, m, imag)); */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   /\* free memory *\/ */
/*   gsl_vector_free(homo_grow_v); */
/*   gsl_vector_free(homo_decay_v); */
/*   gsl_matrix_free(source_matrix); */
/*   gsl_matrix_free(particular_matrix); */
/*   radial_matrix_free(nz, nt, radial_matrix); */
/*   ylm_coeff_free(particular_ylm_coeff); */
/*   ylm_coeff_free(homo_grow_ylm_coeff); */
/*   ylm_coeff_free(homo_decay_ylm_coeff); */
/* } */
  
  
/***************************************************************************************************************/
/* Find the {l, m}th eigenfunction for the radial Poisson equation in terms of Chebyshev coefficients T_i(xi). */
/***************************************************************************************************************/
void solve_radial_particular(int L, gsl_matrix ***radial_matrix, gsl_vector *alpha_vector, gsl_vector *beta_vector, gsl_matrix *source_matrix, gsl_matrix *particular_matrix)
{
  double alpha;
  double beta;
  int z;
  int i;
  int nz;
  int nr;
  gsl_vector *source_vector;
  gsl_vector *particular_vector;
  
  nz = particular_matrix->size1;
  nr = particular_matrix->size2;

  source_vector = gsl_vector_alloc(nr);
  particular_vector = gsl_vector_alloc(nr);

  z = 0;
  /* decompose and solve kernel when L is even: */
  if(!(L%2)){
    alpha = gsl_vector_get(alpha_vector, z);
    gsl_matrix_get_row(source_vector, source_matrix, z);
    solve_kernel_even(L, radial_matrix[z][L], alpha, source_vector, particular_vector);
    gsl_matrix_set_row(particular_matrix, z, particular_vector);
    
    /* decompose and solve kernel when L is odd: */
  } else {
    alpha = gsl_vector_get(alpha_vector, z);
    for(i=0; i<nr; i++) {
      gsl_vector_set(source_vector, i, gsl_matrix_get(source_matrix, z, i));
    }
    solve_kernel_odd(L, radial_matrix[z][L], alpha, source_vector, particular_vector);
    for(i=0; i<nr; i++) {
      gsl_matrix_set(particular_matrix, z, i, gsl_vector_get(particular_vector, i));
    }
  }
  
  /* decompose and solve shell(s) if there are any: */
  for(z=1; z<=nz-2; z++){ /* only evaluated if nz >= 3 */
    alpha = gsl_vector_get(alpha_vector, z);
    beta = gsl_vector_get(beta_vector, z);
    for(i=0; i<nr; i++) {
      gsl_vector_set(source_vector, i, gsl_matrix_get(source_matrix, z, i));
    }
    solve_shell(L, radial_matrix[z][L], alpha, beta, source_vector, particular_vector);
    for(i=0; i<nr; i++) {
      gsl_matrix_set(particular_matrix, z, i, gsl_vector_get(particular_vector, i));
    }
  }
  
  z = nz-1;
  /* decompose and solve external domain: */
  alpha = gsl_vector_get(alpha_vector, z);
  for(i=0; i<nr; i++) {
    gsl_vector_set(source_vector, i, gsl_matrix_get(source_matrix, z, i));
  }
/*   printf("hello\n"); */
/*   print_matrix(source_matrix); */
/*   print_vector(source_vector); */
  solve_ext(L, radial_matrix[z][L], alpha, source_vector, particular_vector);
  for(i=0; i<nr; i++) {
    gsl_matrix_set(particular_matrix, z, i, gsl_vector_get(particular_vector, i));
  }
/*   printf("%f\n", alpha); */
/*   print_matrix(radial_matrix[z][L]); */
/*   print_vector(particular_vector); */
/*   print_matrix(particular_matrix); */
  
  /* free memory */
  gsl_vector_free(source_vector);
  gsl_vector_free(particular_vector);
}


/********************************************************************/
/* Find coefficients for homogeneous solution to satisfy continuity */
/* given the particular solution particular_m.                      */
/* Returns homo_grow_v = (A_0, A_1,..., A_{nz-2}, 0)                */
/* and homo_decay_v = (0, B_1,..., B_{nz-2}, B_{nz-1})              */
/********************************************************************/
void solve_radial_homogeneous(int L, gsl_vector *alpha_v, gsl_vector *beta_v, gsl_matrix *particular_m, gsl_vector *homo_grow_v, gsl_vector *homo_decay_v)
{
  int z;
  int i;
  int nz = particular_m->size1; /* number of zones (must be >= 2) */
  int nr = particular_m->size2;
  
  double R;
  double alpha;

  double fin;
  double dfin;
  double fout;
  double dfout;
  
  gsl_vector *boundary = gsl_vector_calloc (nz-1);
  
  gsl_matrix *cont_m = gsl_matrix_calloc (2*nz-2, 2*nz-2); /* set elements to zero */
  gsl_vector *cont_v = gsl_vector_calloc (2*nz-2);
  gsl_vector *homo_coeff = gsl_vector_calloc (2*nz-2);
  
  int luint;
  gsl_permutation *permute = gsl_permutation_alloc (2*nz-2);
  
  /* nz zones */
  /* nz-1 boundaries between zones */
  /* nz-2 is index of last boundary in boundary vector */
  /* 0 is index of boundary between kernal and 1st shell */
  /* 1...nz-3 are indices of boundaries between two shells */
  /* nz-2 is index of boundary between last shell and external domain */
  
  /* Determine radial positions of the boundaries */
  /* !!This will need to be modified to include the F, G functions in the equation for the boundaries */
  /* when using non spherical boundaries!! */
  for(z=0; z<=nz-2; z++) {
    if(z==0) {
      gsl_vector_set(boundary, 0, gsl_vector_get(alpha_v, 0));
    } else {
      gsl_vector_set(boundary, z, gsl_vector_get(alpha_v, z) + gsl_vector_get(beta_v, z));
    }
  }
  
  
  /*>>>>>>>>>>>>>>> set values for the matrix cont_m in the continuity equation <<<<<<<<<<<<<<*/
  
  if(nz < 2){ /* error */
    
    printf("There must be at least 2 zones");
    
  }else if(nz == 2){ /* there is only a kernel and an external domain */
    
    /* row for continuity */
    gsl_matrix_set(cont_m, 0, 0, pow(gsl_vector_get(boundary, 0), L));
    gsl_matrix_set(cont_m, 0, 1, -pow(gsl_vector_get(boundary, 0), -L-1));
    /* row for continuity of 1st derivative */
    gsl_matrix_set(cont_m, 1, 0, L*pow(gsl_vector_get(boundary, 0), L-1));
    gsl_matrix_set(cont_m, 1, 1, (L+1)*pow(gsl_vector_get(boundary, 0), -L-2));
    
  }else{ /* there are shells */
    
    /* internal boundary: */
    /* row for continuity */
    gsl_matrix_set(cont_m, 0, 0, pow(gsl_vector_get(boundary, 0), L));
    gsl_matrix_set(cont_m, 0, 1, -pow(gsl_vector_get(boundary, 0), L));
    gsl_matrix_set(cont_m, 0, 2, -pow(gsl_vector_get(boundary, 0), -L-1));
    /* row for continuity of 1st derivative */
    gsl_matrix_set(cont_m, 1, 0, L*pow(gsl_vector_get(boundary, 0), L-1));
    gsl_matrix_set(cont_m, 1, 1, -L*pow(gsl_vector_get(boundary, 0), L-1));
    gsl_matrix_set(cont_m, 1, 2, (L+1)*pow(gsl_vector_get(boundary, 0), -L-2));
    
    /* shell boundaries: */
    
    for(z=1; z<=nz-3; z++){
      /* row for continuity */
      gsl_matrix_set(cont_m, 2*z, 2*z-1, pow(gsl_vector_get(boundary, z), L));
      gsl_matrix_set(cont_m, 2*z, 2*z, pow(gsl_vector_get(boundary, z), -L-1)); /* bug: was -L+1 */
      gsl_matrix_set(cont_m, 2*z, 2*z+1, -pow(gsl_vector_get(boundary, z), L));
      gsl_matrix_set(cont_m, 2*z, 2*z+2, -pow(gsl_vector_get(boundary, z), -L-1));
      /* row for continuity of 1st derivative */
      gsl_matrix_set(cont_m, 2*z+1, 2*z-1, L*pow(gsl_vector_get(boundary, z), L-1));
      gsl_matrix_set(cont_m, 2*z+1, 2*z, -(L+1)*pow(gsl_vector_get(boundary, z), -L-2));
      gsl_matrix_set(cont_m, 2*z+1, 2*z+1, -L*pow(gsl_vector_get(boundary, z), L-1));
      gsl_matrix_set(cont_m, 2*z+1, 2*z+2, (L+1)*pow(gsl_vector_get(boundary, z), -L-2));
    }
    
    /* external boundary: */
    /* row for continuity */
    gsl_matrix_set(cont_m, 2*nz-4, 2*nz-5, pow(gsl_vector_get(boundary, nz-2), L));
    gsl_matrix_set(cont_m, 2*nz-4, 2*nz-4, pow(gsl_vector_get(boundary, nz-2), -L-1));
    gsl_matrix_set(cont_m, 2*nz-4, 2*nz-3, -pow(gsl_vector_get(boundary, nz-2), -L-1));
    /* row for continuity of 1st derivative */
    gsl_matrix_set(cont_m, 2*nz-3, 2*nz-5, L*pow(gsl_vector_get(boundary, nz-2), L-1));
    gsl_matrix_set(cont_m, 2*nz-3, 2*nz-4, -(L+1)*pow(gsl_vector_get(boundary, nz-2), -L-2));
    gsl_matrix_set(cont_m, 2*nz-3, 2*nz-3, (L+1)*pow(gsl_vector_get(boundary, nz-2), -L-2));
    
  }
  
  /*print_matrix(cont_m);*/
  
  /*>>>>>>>>>>>>>>> set values for the vector cont_v in the continuity equation <<<<<<<<<<<<<<*/

  /*!!!!!THE FOLLOWING DOES NOT SUM FROM SMALLEST TO GREATEST!!!!!! It could matter*/

  if(nz == 2){ /* there is only a kernel and an external domain */
    
    fin = dfin = fout = dfout = 0;
    R = gsl_vector_get(boundary, 0);
    /* function and derivative on inside: */
    alpha = R;
    if(!(L%2)){ /* L is even */
      for(i=0; i<nr; i++){
	fin += gsl_matrix_get(particular_m, 0, i);
	dfin += 4*i*i*gsl_matrix_get(particular_m, 0, i);
      }
      dfin /= alpha;
    } else { /* L is odd */
      for(i=0; i<nr; i++){
	fin += gsl_matrix_get(particular_m, 0, i);
	dfin += (2*i+1)*(2*i+1)*gsl_matrix_get(particular_m, 0, i);
      }
      dfin /= alpha;
    }
    /* function and derivative on outside: */
    alpha = -0.5/R;
    for(i=0; i<nr; i++){
      fout += neg1toi(i)*gsl_matrix_get(particular_m, 1, i);
      dfout += neg1toi(i+1)*i*i*gsl_matrix_get(particular_m, 1, i);
    }
    dfout /= -(alpha*R*R);
    
    gsl_vector_set(cont_v, 0, -fin+fout);
    gsl_vector_set(cont_v, 1, -dfin+dfout);
    
  }else{ /* there are shells */
    
    /* interface between kernel and first shell */
    fin = dfin = fout = dfout = 0;
    /* function and derivative on inside: */
    R = gsl_vector_get(boundary, 0);
    alpha = R;
    if(!(L%2)){ /* L is even */
      for(i=0; i<nr; i++){
	fin += gsl_matrix_get(particular_m, 0, i);
	dfin += 4*i*i*gsl_matrix_get(particular_m, 0, i);
      }
      dfin /= alpha;
    } else { /* L is odd */
      for(i=0; i<nr; i++){
	fin += gsl_matrix_get(particular_m, 0, i);
	dfin += (2*i+1)*(2*i+1)*gsl_matrix_get(particular_m, 0, i);
      }
      dfin /= alpha;
    }
    /* function and derivative on outside: */
    alpha = 0.5*(gsl_vector_get(boundary, 1)-R);
    for(i=0; i<nr; i++){
      fout += neg1toi(i)*gsl_matrix_get(particular_m, 1, i);
      dfout += neg1toi(i+1)*i*i*gsl_matrix_get(particular_m, 1, i);
    }
    dfout /= alpha;
    
    gsl_vector_set(cont_v, 0, -fin+fout);
    gsl_vector_set(cont_v, 1, -dfin+dfout);
    
    /* interface between shells */
    /* (skipped if nz<=3) */
    for(z=1; z<=nz-3; z++){
      fin = dfin = fout = dfout = 0;
      R = gsl_vector_get(boundary, z);
      /* function and derivative on inside: */
      alpha = 0.5*(R-gsl_vector_get(boundary, z-1));
      for(i=0; i<nr; i++){
	fin += gsl_matrix_get(particular_m, z, i);
	dfin += i*i*gsl_matrix_get(particular_m, z, i);
      }
      dfin /= alpha;
      /* function and derivative on outside: */
      alpha = 0.5*(gsl_vector_get(boundary, z+1)-R);
      for(i=0; i<nr; i++){
	fout += neg1toi(i)*gsl_matrix_get(particular_m, z+1, i);
	dfout += neg1toi(i+1)*i*i*gsl_matrix_get(particular_m, z+1, i);
      }
      dfout /= alpha;
      
      gsl_vector_set(cont_v, 2*z, -fin+fout);
      gsl_vector_set(cont_v, 2*z+1, -dfin+dfout);
    }
    
    /* interface between last shell and external domain */
    fin = dfin = fout = dfout = 0;
    R = gsl_vector_get(boundary, nz-2);
    /* function and derivative on inside: */
    alpha = 0.5*(R-gsl_vector_get(boundary, nz-3));
    for(i=0; i<nr; i++){
      fin += gsl_matrix_get(particular_m, nz-2, i);
      dfin += i*i*gsl_matrix_get(particular_m, nz-2, i);
    }
    dfin /= alpha;
    /* function and derivative on outside: */
    alpha = -0.5/R;
    for(i=0; i<nr; i++){
      fout += neg1toi(i)*gsl_matrix_get(particular_m, nz-1, i);
      dfout += neg1toi(i+1)*i*i*gsl_matrix_get(particular_m, nz-1, i);
    }
    dfout /= -(alpha*R*R);
    
    gsl_vector_set(cont_v, 2*nz-4, -fin+fout);
    gsl_vector_set(cont_v, 2*nz-3, -dfin+dfout);
  }
  
  /*print_vector(cont_v);*/
  
  /*>>>>>>>>>>>>>> Solve equation for continuity <<<<<<<<<<<<<<<<<*/
  
  /* Solve for homo_coeff in cont_m.homo_coeff = cont_v */
  gsl_linalg_LU_decomp (cont_m, permute, &luint);
  gsl_linalg_LU_solve (cont_m, permute, cont_v, homo_coeff);
  
  /* Place elements of homo_coeff into homo_grow_v, homo_decay_v */
  
  gsl_vector_set(homo_grow_v, 0, gsl_vector_get(homo_coeff, 0));
  gsl_vector_set(homo_decay_v, 0, 0.0);
  for(z=1; z<=nz-2; z++) {
    gsl_vector_set(homo_grow_v, z, gsl_vector_get(homo_coeff, 2*z-1));
    gsl_vector_set(homo_decay_v, z, gsl_vector_get(homo_coeff, 2*z));
  }
  gsl_vector_set(homo_grow_v, nz-1, 0.0);
  gsl_vector_set(homo_decay_v, nz-1, gsl_vector_get(homo_coeff, 2*nz-3));
  
  /* free memory */
  gsl_vector_free(boundary);
  gsl_matrix_free(cont_m);
  gsl_vector_free(cont_v);
  gsl_vector_free(homo_coeff);
  gsl_permutation_free(permute);
}


/*******************************************************************************/
/* Take coefficients {A_Lm, B_Lm} of homogeneous solution for each {L, m} and  */
/* return the coefficients C_iLm of T_i(xi).                                   */
/* f_h(r) = A_Lm r^L + B_Lm r^-(L+1)  -->  f_h(xi) = C_iLm T_i(xi)             */
/*******************************************************************************/
void homogeneoustochebyshev(ylm_coeff *homo_grow_ylm_coeff, ylm_coeff *homo_decay_ylm_coeff, gsl_vector *alpha_vector, gsl_vector *beta_vector, ylm_coeff *homogeneous_ylm_coeff)
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction */
  int nt; /* number of points in theta direction */
  int z, i, L, m;
  int imag;
  int ireverse;
  double alpha, beta;
  double alm, blm;
  double xi_i;
  double *in1dxi; /* even part of kernel or all parts of other zones */
  double *in1dxiodd; /* odd part of kernel */
  fftw_plan plan_forward_xi;
  fftw_plan plan_forward_xiodd;
  double *out1dxi;
  double *out1dxiodd;
  
  nz = homogeneous_ylm_coeff->nz;
  nr = homogeneous_ylm_coeff->nr;
  nt = homogeneous_ylm_coeff->nt;
    
  /* Set up arrays to hold the data */
  in1dxi = fftw_malloc ( sizeof ( double ) * nr );
  in1dxiodd = fftw_malloc ( sizeof ( double ) * (nr-1) );
  out1dxi = fftw_malloc ( sizeof ( double ) * nr );
  out1dxiodd = fftw_malloc ( sizeof ( double ) * (nr-1) );
  
  /* Create FFTW plans.  The algorithm will be stored in the fftw_plan structure.         */
  /* The arrays can be changed later (must stay same size) but the plan will be the same. */
  /* The plan usually overwrites the data, so set the data after making the plan.         */
  /* It might be faster to figure out how to use FFTW wisdom to speed up planning.        */
  plan_forward_xi = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
  /* do type-3 DCT (inverse of type-2 DCT) */ 
  plan_forward_xiodd = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT01, FFTW_ESTIMATE ); 
 
/*   z = 0; /\* kernel *\/ */
/*   /\* even Chebyshev series if L is even *\/  */
/*   for ( imag=0; imag<=1; imag++ ) { */
/*     for ( m = imag; m < nt-imag; m++ ) { */
/*       for ( L = m; L < nt-m%2; L++ ) { */
/* 	alpha = gsl_vector_get(alpha_vector, z); */
/* 	alm = ylm_coeff_get(homo_grow_ylm_coeff, z, 0, L, m, imag);   */
/* 	for ( i = 0; i < nr; i++ ) { */
/* 	  xi_i = sin(PI*i/(2*(nr-1))); */
/* 	  in1dxi[i] = alm*pow(alpha*xi_i, L); */
/* 	} */
/* 	fftw_execute ( plan_forward_xi ); */
/* 	for ( i = 0; i < nr; i++ ) { */
/* 	  ylm_coeff_set( homogeneous_ylm_coeff, z, i, L, m, imag,  */
/* 			 out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) ); */
/* 	} */
/*       } */
/*     } */
/*   } */

  z = 0; /* kernel */
  /* even Chebyshev series if L is even */
  for ( imag=0; imag<=1; imag++ ) {
    for ( L = 2*imag ; L < nt-imag; L += 2 ) { /* L = 0 and L = npc-1 don't have immaginary parts */
      for ( m = 0; m <= L; m++ ) {
	alpha = gsl_vector_get(alpha_vector, z);
	alm = ylm_coeff_get(homo_grow_ylm_coeff, z, 0, L, m, imag);
	for ( i = 0; i < nr; i++ ) {
	  xi_i = sin(PI*i/(2*(nr-1)));
	  in1dxi[i] = alm*pow(alpha*xi_i, L);
	}
	fftw_execute ( plan_forward_xi );
	for ( i = 0; i < nr; i++ ) {
	  ylm_coeff_set( homogeneous_ylm_coeff, z, i, L, m, imag,
		     out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
	}
      }
    }
  }
    
  z = 0; /* kernel */
  /* odd Chebyshev series if L is odd */
  for ( imag=0; imag<=1; imag++ ) {
    for ( L = 1 ; L < nt-imag; L += 2 ) { /* L = npc-1 doesn't have immaginary part */
      for ( m = 0; m <= L; m++ ) {    
	alpha = gsl_vector_get(alpha_vector, z);
	alm = ylm_coeff_get(homo_grow_ylm_coeff, z, 0, L, m, imag);  
	/*printf("in: imag=%d, L=%d, m=%d\n", imag, L, m);*/
	for ( i = 1; i < nr; i++ ) {
	  ireverse = nr - i;
	  xi_i = sin(PI*ireverse/(2*(nr-1)));
	  in1dxiodd[i-1] = alm*pow(alpha*xi_i, L); 
	  /*printf("%f\t", in1dxiodd[i]);*/
	}
	/*printf("\nout:");*/
	/* do type-3 DCT */
	fftw_execute ( plan_forward_xiodd );
	for ( i = 0; i < nr-1; i++ ) {
	  ylm_coeff_set( homogeneous_ylm_coeff, z, i, L, m, imag, 
			 out1dxiodd[i]/(nr-1) );
	  /*printf("%f\t", out1dxiodd[i]);*/
	}	
	/*printf("\n");*/
	ylm_coeff_set( homogeneous_ylm_coeff, z, nr-1, L, m, imag, 0.0 ); /* set coefficient for T_{nr-1}(xi) = 0 */
      }
    }
  }
      
  /* other zones */
  for(z=1; z<=nz-2; z++){ /* only evaluated if nz >= 3 */
    /* normal Chebyshev series */
    for ( imag=0; imag<=1; imag++ ) {
      for ( L = imag ; L < nt-imag; L++ ) { /* L = 0 and L = npc-1 don't have immaginary parts */
	for ( m = 0; m <= L; m++ ) { 
	  alpha = gsl_vector_get(alpha_vector, z);
	  beta = gsl_vector_get(beta_vector, z);
	  alm = ylm_coeff_get(homo_grow_ylm_coeff, z, 0, L, m, imag);
	  blm = ylm_coeff_get(homo_decay_ylm_coeff, z, 0, L, m, imag);     
	  for ( i = 0; i < nr; i++ ) {
	    xi_i = -cos(PI*i/(nr-1));
	    in1dxi[i] = alm*pow(alpha*xi_i + beta, L) + blm*pow(alpha*xi_i + beta, -L-1);
	  }
	  fftw_execute ( plan_forward_xi );
	  for ( i = 0; i < nr; i++ ) {
	    ylm_coeff_set( homogeneous_ylm_coeff, z, i, L, m, imag, 
			   out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
	  }
	}
      }
    }
  }
  
  /* external zone */
  z = nz-1;  
  /* normal Chebyshev series */
  for ( imag=0; imag<=1; imag++ ) {
    for ( L = imag ; L < nt-imag; L++ ) { /* L = 0 and L = npc-1 don't have immaginary parts */
      for ( m = 0; m <= L; m++ ) { 
	alpha = gsl_vector_get(alpha_vector, z);
	blm = ylm_coeff_get(homo_decay_ylm_coeff, z, 0, L, m, imag);     
	for ( i = 0; i < nr; i++ ) {
	  xi_i = -cos(PI*i/(nr-1));
	  in1dxi[i] = blm*pow(alpha*(xi_i-1), L+1);
	}
	fftw_execute ( plan_forward_xi );
	for ( i = 0; i < nr; i++ ) {
	  ylm_coeff_set( homogeneous_ylm_coeff, z, i, L, m, imag, 
			 out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
	}
      }
    }
  }
  
  /* Delete plan and arrays */
  fftw_free ( in1dxi );
  fftw_free ( out1dxi );
  fftw_free ( in1dxiodd );
  fftw_free ( out1dxiodd );
}



/* /\*******************************************************************\/ */
/* /\* Take coefficients a, b of the homogeneous part of the solution  *\/ */
/* /\* f_homo = a r^L + b r^{-L-1}                                     *\/ */
/* /\* and find the corresponding Chebyshev coefficients T_i(xi)       *\/ */
/* /\*******************************************************************\/ */
/* void homo_to_cheb_coeff_kernel_even(int L, double a, double alpha, gsl_vector *cheb_coeff_v); */

/* void homo_to_cheb_coeff_kernel_odd(int L, double a, double alpha, gsl_vector *cheb_coeff_v); */

/* void homo_to_cheb_coeff_shell(int L, double a, double b, double alpha, double beta, gsl_vector *cheb_coeff_v) */
/* { */
/*   int i; */
/*   int nr; */
/*   double *in1dxi; */
/*   fftw_plan plan_forward_xi; */
/*   double *out1dxi; */
  
/*   nr = cheb_coeff_v->size; */

/*   in1dxi = fftw_malloc ( sizeof ( double ) * nr ); */
/*   out1dxi = fftw_malloc ( sizeof ( double ) * nr ); */
/*   plan_forward_xi = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE ); */

/*   for ( i = 0; i < nr; i++ ) { */
/*     xi = -cos(PI*i/(nr-1)); */
/*     r = alpha*(xi  */
/* 	       + 0.25*(xi*xi*xi-3.0*xi+2.0)*f_lm */
/* 	       + 0.25*(-xi*xi*xi+3.0*xi+2.0)*g_lm) */
/*       + beta; */
/*     in1dxi[i] = a*pow(r, L) + b*pow(r, -L-1); */
/*   } */
/*   fftw_execute ( plan_forward_xi ); */
/*   for ( i = 0; i < nr; i++ ) { */
/*     gsl_vector_set(cheb_coeff_v, i, out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) ); */
/*   } */
  
/*  fftw_destroy_plan ( plan_forward_xi ); */
/* } */

/* void homo_to_cheb_coeff_ext(int L, double b, double alpha, gsl_vector *cheb_coeff_v); */



/*******************************************************/
/* Solve even part of kernel:                          */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_kernel_even(int L, gsl_matrix *A_even_source, double alpha, gsl_vector *s_even, gsl_vector *f_even)
{
  int N = (*s_even).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */
  
  gsl_matrix *A_even;
  gsl_matrix *A_even_trunc;

  gsl_vector *s_even_trunc;
  gsl_vector *f_even_trunc;

  int k;
  gsl_permutation *permute;

  double fzero; /* f(x=0) */

  
  if(L%2) {
	printf("L is not even in solve_kernel_even.");
	p = 2;
  } else if(L == 0)
	p = 1;
  else
	p = 2;
  
  
  gsl_vector_set_zero(f_even);
  
  A_even = gsl_matrix_alloc (N, N);
  A_even_trunc = gsl_matrix_alloc (N-p, N-p);
  
  s_even_trunc = gsl_vector_alloc (N-p);
  f_even_trunc = gsl_vector_alloc (N-p);
  
  permute = gsl_permutation_alloc (N-p);
  
  /* copy operator matrix so original is not destroyed */
  gsl_matrix_memcpy(A_even, A_even_source);

  /* make operator matrix */
  /*set_A_kernel_even(A_even, L);*/
  /*printf("A_even:\n");
  print_matrix(A_even);
  print_vector(s_even);*/

  /* scale s_kernel_even by \alpha^2 */
  gsl_vector_scale(s_even, alpha*alpha);
  
  /* make matrix banded and perform the same operations on source vector */
  /* L_i = (1+\delta_{0i})L_i - L_{i+2}, for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_even, i, j,
					 (1+delta(i, 0))*gsl_matrix_get(A_even, i, j) - gsl_matrix_get(A_even, i+2, j));
	gsl_vector_set(s_even, i,
				   (1+delta(i, 0))*gsl_vector_get(s_even, i) - gsl_vector_get(s_even, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_even, i, j, gsl_matrix_get(A_even, i, j) - gsl_matrix_get(A_even, i+2, j));
	gsl_vector_set(s_even, i,
				   gsl_vector_get(s_even, i) - gsl_vector_get(s_even, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_even, i, j, gsl_matrix_get(A_even, i, j) - gsl_matrix_get(A_even, i+1, j));
	gsl_vector_set(s_even, i,
				   gsl_vector_get(s_even, i) - gsl_vector_get(s_even, i+1));
  }
  /*printf("A_even banded:\n");
  print_matrix(A_even);
  print_vector(s_even);*/
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of source vector                      */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_even_trunc, i, j, gsl_matrix_get(A_even, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_even_trunc, i, gsl_vector_get(s_even, i));
  /*printf("A_even banded truncated:\n");
  print_matrix(A_even_trunc);
  print_vector(s_even_trunc);*/
  
  /********************************************************/
  /* TO DO: Use LAPACK to find f instead of what's below, */
  /*        to take advantage of the banded matrix        */
  /********************************************************/

  /* Solve for f_even in A_even_trunc.f_even_trunc = s_even_trunc */
  gsl_linalg_LU_decomp (A_even_trunc, permute, &k);
  gsl_linalg_LU_solve (A_even_trunc, permute, s_even_trunc, f_even_trunc);
  /*print_vector(f_even_trunc);*/
  /*printf("f_even_trunc:\n");
	print_vector(f_even_trunc);*/


  /* Set solution vector f_even = (0, ..., f_even_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_even, i+p, gsl_vector_get(f_even_trunc, i));
  /*printf("f_even after restoring leading zeros:\n");
	print_vector(f_even);*/

  /* satisfy boundary conditions */
  /* f(x=r=0)=0 for L>=2 */
  /* df/dr(r=0)=0 is already satisfied for even L */
  if(L>=2){
	fzero = 0.0;
	for(i=0; i<N; i++)
	  fzero += neg1toi(i)*gsl_vector_get(f_even, i);
	gsl_vector_set(f_even, 0, -fzero); /* first element of f was already zero */
  }
  /*printf("f_even after requiring solution equal 0 at x=0:\n");
	print_vector(f_even);*/

  /* free memory */
  gsl_matrix_free(A_even);
  gsl_matrix_free(A_even_trunc);
  gsl_vector_free(s_even_trunc);
  gsl_vector_free(f_even_trunc);
  gsl_permutation_free(permute);
}


/*******************************************************/
/* Solve odd part of kernel:                           */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_kernel_odd(int L, gsl_matrix *A_odd_source, double alpha, gsl_vector *s_odd, gsl_vector *f_odd)
{
  int N = (*s_odd).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */
  
  gsl_matrix *A_odd;
  gsl_matrix *A_odd_trunc;
  
  gsl_vector *s_odd_trunc;
  gsl_vector *f_odd_trunc;
  
  int k;
  gsl_permutation *permute;
  
  double dfzero;
  
  if(!(L%2)) {
	printf("L is not odd in solve_kernel_odd.");
	p = 2;
  } else if(L == 1)
	p = 1;
  else
	p = 2;
  
  gsl_vector_set_zero(f_odd);
  A_odd = gsl_matrix_alloc (N, N);
  A_odd_trunc = gsl_matrix_alloc (N-p, N-p);
  
  s_odd_trunc = gsl_vector_alloc (N-p);
  f_odd_trunc = gsl_vector_alloc (N-p);
  
  permute = gsl_permutation_alloc (N-p);
  
  /* copy operator matrix so original is not destroyed */
  gsl_matrix_memcpy(A_odd, A_odd_source);

  /* make operator matrix */
  /*set_A_kernel_odd(A_odd, L);*/
  /*printf("A_kernel_odd:\n");
	print_matrix(A_odd);	*/
  
  /* scale s_kernel_even by alpha^2 */
  gsl_vector_scale(s_odd, alpha*alpha);
  /*print_vector(s_odd);*/

  /* make matrix banded and perform the same operations on source vector */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j,gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+2, j));
	gsl_vector_set(s_odd, i,
				   gsl_vector_get(s_odd, i) - gsl_vector_get(s_odd, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+2, j));
	gsl_vector_set(s_odd, i,
				   gsl_vector_get(s_odd, i) - gsl_vector_get(s_odd, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+1, j));
	gsl_vector_set(s_odd, i,
				   gsl_vector_get(s_odd, i) - gsl_vector_get(s_odd, i+1));
  }
  /*printf("A_kernel_odd banded:\n");
  print_matrix(A_odd);
  print_vector(s_odd);*/
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_odd_trunc, i, j, gsl_matrix_get(A_odd, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_odd_trunc, i, gsl_vector_get(s_odd, i));
   
  /*printf("A_kernel_odd banded truncated:\n");
  print_matrix(A_odd_trunc);
  print_vector(s_odd_trunc);*/
  
  /* Solve for f_odd in A_odd_trunc.f_odd_trunc = s_odd_trunc */
  gsl_linalg_LU_decomp (A_odd_trunc, permute, &k);
  gsl_linalg_LU_solve (A_odd_trunc, permute, s_odd_trunc, f_odd_trunc);
  /*printf("f_odd_trunc:\n");
	print_vector(f_odd_trunc);*/

  /* Set solution vector f_odd = (0, ..., f_odd_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_odd, i+p, gsl_vector_get(f_odd_trunc, i));
  /*printf("f_odd after restoring leading zeros:\n");
	print_vector(f_odd);*/

  /* satisfy boundary conditions */
  /* df/dr(r=0)=0 for L>=3 */
  /* f(x=r=0)=0 is already satisfied for odd L */
  if(L>=3){
	dfzero = 0.0;
	for(i=0; i<N; i++){
	  dfzero += neg1toi(i)*(2*i+1)*gsl_vector_get(f_odd, i);
	  /*printf("%f\n", dfzero);*/
	}
	gsl_vector_set(f_odd, 0, -dfzero); /* first element of f was already zero */
  }
  /*printf("f_odd after making slope 0 at x=0:\n");
	print_vector(f_odd);*/
  
  /* free memory */
  gsl_matrix_free(A_odd);
  gsl_matrix_free(A_odd_trunc);
  gsl_vector_free(s_odd_trunc);
  gsl_vector_free(f_odd_trunc);
  gsl_permutation_free(permute);

}


/*******************************************************/
/* Solve shell:                                        */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_shell(int L, gsl_matrix *A_shell_source, double alpha, double beta, gsl_vector *s_old, gsl_vector *f)
{
  int N = (*s_old).size;
  int i, j;
  int p=2; /* order of eigenvalue for homogeneous solution to A */
  
  gsl_matrix *A = gsl_matrix_alloc (N, N);
  gsl_matrix *x = gsl_matrix_alloc (N, N);
  gsl_matrix *id = gsl_matrix_alloc (N, N);
  gsl_matrix *axplusb = gsl_matrix_calloc (N, N); /* set to zero */
  gsl_matrix *axplusb2 = gsl_matrix_alloc (N, N);

  gsl_vector *s = gsl_vector_alloc (N);
 
  gsl_matrix *A_trunc = gsl_matrix_alloc (N-p, N-p);
 
  gsl_vector *hello = gsl_vector_calloc (N-p);

  gsl_vector *s_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_trunc = gsl_vector_alloc (N-p);
  
  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);
  
  gsl_vector_set_zero(f);

  /*printf("alpha=%f\tbeta=%f\n", alpha, beta);*/
  
  /* copy operator matrix so original is not destroyed */
  gsl_matrix_memcpy(A, A_shell_source);

  /* make operator matrix */
  /*set_A_shell(A, L, alpha, beta);*/
  /*set_A_shell(A, 7, 3.345, 73.2345);*/
  /*gsl_matrix_set_identity(A);*/
  /*printf("A shell:\n");
	print_matrix(A);*/
  
  /*printf("original source before multiplying by (alpha x + beta)^2:\n");
	print_vector(s_old);*/
  
  /* multiply s by (alpha x + beta)^2 */
  set_x(x);
  gsl_matrix_scale(x, alpha);
  gsl_matrix_set_identity(id);
  gsl_matrix_scale(id, beta);
  gsl_matrix_add(axplusb, x);
  gsl_matrix_add(axplusb, id);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, axplusb, axplusb,
				 0.0, axplusb2);
  gsl_blas_dgemv(CblasNoTrans,
				 1.0, axplusb2, s_old,
				 0.0, s);
  /*printf("source after multiplying by (alpha x + beta)^2:\n");
	print_vector(s);*/

  /* make matrix banded and perform the same operations on source vector */
  /* L_i = \frac{(1+\delta_{0i})L_i - L_{i+2}}{i+1}, for 0 <= i <= N-3   */
  /* L_i = L_i - L_{i+2},                            for 0 <= i <= N-5   */
  /* for(i=0; i<=N-3; i++) { */
/* 	for(j=0; j<N; j++) */
/* 	  gsl_matrix_set(A, i, j, */
/* 					 ((1+delta(i, 0))*gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j))/(i+1)); */
/* 	gsl_vector_set(s, i, */
/* 				   ((1+delta(i, 0))*gsl_vector_get(s, i) - gsl_vector_get(s, i+2))/(i+1)); */
/*   } */
/*   for(i=0; i<=N-5; i++) { */
/* 	for(j=0; j<N; j++) */
/* 	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j)); */
/* 	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2)); */
/*   } */
/*   printf("A shell banded:\n"); */
/*   print_matrix(A); */
/*   print_vector(s); */
  

  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_trunc, i, j, gsl_matrix_get(A, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_trunc, i, gsl_vector_get(s, i));
  

  /*gsl_matrix_set_identity(A_trunc);*/


  /*printf("A shell banded truncated:\n");
  print_matrix(A_trunc);
  print_vector(s_trunc);*/
  
  
  
  /* Solve for f_trunc in A_trunc.f_trunc = s_trunc */
  gsl_linalg_LU_decomp (A_trunc, permute, &k);
  gsl_linalg_LU_solve (A_trunc, permute, s_trunc, f_trunc);
  /*printf("f_trunc for shell:\n");
	print_vector(f_trunc);*/
  

  
  gsl_blas_dgemv(CblasNoTrans,
				 1.0, A_trunc, f_trunc,
				 0.0, hello);
  /*printf("A_trunc.f_trunc should equal s_trunc:\n");
	print_vector(hello);*/
  



  /* Set solution vector f = (0, ..., f_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f, i+p, gsl_vector_get(f_trunc, i));
  /*printf("f after restoring leading zeros:\n");
	print_vector(f);*/
  

  /* free memory */
  gsl_matrix_free(A);
  gsl_matrix_free(x);
  gsl_matrix_free(id);
  gsl_matrix_free(axplusb);
  gsl_matrix_free(axplusb2);
  gsl_matrix_free(A_trunc);
  
  gsl_vector_free(s);
  gsl_vector_free(s_trunc);
  gsl_vector_free(f_trunc);
    
  gsl_permutation_free(permute);
}


/*******************************************************/
/* Solve external domain:                              */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_ext(int L, gsl_matrix *A_ext_source, double alpha, gsl_vector *s_old, gsl_vector *f)
{
  int N = (*s_old).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */
  
  gsl_matrix *A;
  gsl_matrix *xmin1inv;
  gsl_matrix *xmin1inv2;
  gsl_matrix *xmin1inv4;
  gsl_matrix *A_trunc;
  gsl_vector *s;
  gsl_vector *s_trunc;
  gsl_vector *f_trunc;
  
  int k;
  gsl_permutation *permute;
  
  double finf;
  double dfinf;

  if(L == 0)
	p = 2;
  else
	p = 3;
  
  gsl_vector_set_zero(f);
  
  A = gsl_matrix_alloc (N, N);
  xmin1inv = gsl_matrix_alloc (N, N);
  xmin1inv2 = gsl_matrix_alloc (N, N);
  xmin1inv4 = gsl_matrix_alloc (N, N);
  A_trunc = gsl_matrix_alloc (N-p, N-p);
  s = gsl_vector_alloc (N);
  s_trunc = gsl_vector_alloc (N-p);
  f_trunc = gsl_vector_alloc (N-p);
  permute = gsl_permutation_alloc (N-p);
  
  /* copy operator matrix so original is not destroyed */
  gsl_matrix_memcpy(A, A_ext_source);

  /* make operator matrix */
  /*set_A_ext(A, L);*/
  /*printf("A_ext:\n");
	print_matrix(A);*/
  
  /* divide s by \alpha^2(x-1)^4 */
  set_xmin1inv(xmin1inv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, xmin1inv, xmin1inv,
				 0.0, xmin1inv2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, xmin1inv2, xmin1inv2,
				 0.0, xmin1inv4);
  gsl_blas_dgemv(CblasNoTrans,
				 1.0, xmin1inv4, s_old,
				 0.0, s);
  gsl_vector_scale(s, 1/(alpha*alpha));
  
  /*printf("LHS of equation for external domain:\n");
	print_vector(s);*/

  /* make matrix banded and perform the same operations on source vector */
  /* L_i = (1+\delta_{0i})L_i - L_{i+2}, for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j,
					 (1+delta(i, 0))*gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, (1+delta(i, 0))*gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+1, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+1));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
  }
  /*printf("A banded:\n");
  print_matrix(A);
  print_vector(s);*/
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_trunc, i, j, gsl_matrix_get(A, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_trunc, i, gsl_vector_get(s, i));

  /*printf("A banded truncated:\n");
  print_matrix(A_trunc);
  print_vector(s_trunc);*/
  
  /* Solve for f_even in A_even_trunc.f_even_trunc = s_even_trunc */
  gsl_linalg_LU_decomp (A_trunc, permute, &k);
  gsl_linalg_LU_solve (A_trunc, permute, s_trunc, f_trunc);
  /*print_vector(f_trunc);*/
  
  /* Set solution vector f = (0, ..., f_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f, i+p, gsl_vector_get(f_trunc, i));
  /*print_vector(f);*/
 
  /* satisfy boundary conditions */
  /* f_particular(x=1) + \alpha*T_0(x=1) = 0 */
  finf = 0;
  for(i=2; i<N; i++)
	finf += gsl_vector_get(f, i);
  gsl_vector_set(f, 0, -finf); /* first element of f was already zero */
  /* set df/dx=0 for L>=1 by changing T_0 term */
  if(L>=1){
	dfinf = 0;
	for(i=2; i<N; i++)
	  dfinf += i*i*gsl_vector_get(f, i);
	gsl_vector_set(f, 1, -dfinf);
  }
    
  /*printf("external solution before homogeneous part is added\n");
  printf("but after T_0, T_1 terms are accounted for:\n");
  print_vector(f);*/
  
  /* free memory */
  gsl_matrix_free(A);
  gsl_matrix_free(xmin1inv);
  gsl_matrix_free(xmin1inv2);
  gsl_matrix_free(xmin1inv4);
  gsl_matrix_free(A_trunc);
  gsl_vector_free(s);
  gsl_vector_free(s_trunc);
  gsl_vector_free(f_trunc);
  gsl_permutation_free(permute);
}




/************************************************************************************/
/* takes derivative df/dx when f is expressed as coefficients of a Chebyshev series */
/* type represents regular, even, or odd series                                     */
/* f_vector is the coefficients                                                     */
/* xi must be -1 or +1                                                              */
/************************************************************************************/
/*  double dfdxi_edge(int type, gsl_vector *f_vector, int xi) */
/*  { */
   
/*  } */
