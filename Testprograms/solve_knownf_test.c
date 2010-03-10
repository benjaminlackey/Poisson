/*********************************************************************************
 *      *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/residual.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/poisson.h solve_knownf_test.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>

/* fftw header */
#include <fftw3.h>

/* own header */
#include "poisson.h"

void pseudo_spherical_laplacian(scalar3d *pseudolapf_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector);
void d2fdxi2(coeff *d2fdxi2, coeff *f);
double boundary(int z, double theta, double phi);
double source(int z, double r, double theta, double phi);
double field(int z, double r, double theta, double phi);

int main (void)
{
  int z, i, j, k, imag;
  int L, m;
  int nz = 4;
  int nr = 35; /* must be odd? */
  int nt;
  int np = 12; /* must be even */

  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *field_scalar3d;
  scalar3d *source_scalar3d;
  scalar3d *roru_scalar3d;

  scalar3d *residual_scalar3d;
  scalar3d *j1_scalar3d;
  scalar3d *j2_scalar3d;
  scalar3d *j3_scalar3d;
  scalar3d *source_eff_scalar3d;
  scalar3d *pseudolapf_scalar3d;
  scalar3d *field_solution_scalar3d;
  coeff *field_solution_coeff;

  bound_coeff *f_bound_coeff;
  bound_coeff *g_bound_coeff;
  coeff *source_coeff;
  coeff *field_coeff;

  coeff *residual_coeff;
  ylm_coeff *residual_ylm_coeff;
  scalar3d *a_scalar3d;
  coeff *a_coeff;
  ylm_coeff *a_ylm_coeff;

  ylm_coeff *field_ylm_coeff;
  ylm_coeff *field_solution_ylm_coeff;

  gsl_matrix **fouriertoylm_matrix;
  gsl_matrix **ylmtofourier_matrix;

  double roru, xi_i, theta_j, phi_k;
  double alpha, j1, j2, j3, a, residual_d, source_d, source_eff;
  double pseudolapf;
  double num, anal, error;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  source_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  roru_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  residual_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j3_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  source_eff_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  pseudolapf_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_solution_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_solution_coeff = coeff_alloc(nz, nr, nt, np);

  f_bound_coeff = bound_coeff_alloc(nz, nt, np);
  g_bound_coeff = bound_coeff_alloc(nz, nt, np);
  source_coeff = coeff_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);

  residual_coeff = coeff_alloc(nz, nr, nt, np);
  residual_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  a_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  a_coeff = coeff_alloc(nz, nr, nt, np);
  a_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);

  field_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  field_solution_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);

  /* allocate and make matrices for fourier <--> spherical harmonic transforms */
  fouriertoylm_matrix = fouriertoylm_matrix_alloc(nt);
  ylmtofourier_matrix = ylmtofourier_matrix_alloc(nt);
  fouriertoylm_matrix_set(fouriertoylm_matrix);
  ylmtofourier_matrix_set(ylmtofourier_matrix);
  
  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  /* print_scalar2d(boundary_scalar2d); */

  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  print_vector(alpha_vector); 
  print_vector(beta_vector); 
  gridtofourier_bound(f_bound_coeff, f_scalar2d);
  gridtofourier_bound(g_bound_coeff, g_scalar2d);
/*   print_bound_coeff(f_bound_coeff); */
/*   print_bound_coeff(g_bound_coeff); */

  /* evaluate source at gridpoints */
  functiontogrid(source_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, source);
  gridtofourier(source_coeff, source_scalar3d, 0, 0);
/*   print_coeff(source_coeff); */
  
  /* evaluate analytical solution at gridpoints */
  functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, field);
  gridtofourier(field_coeff, field_scalar3d, 0, 0);
/*   print_coeff(field_coeff); */
  
  /* evaluate radius at gridpoints */
  rofxtp(roru_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* evaluate residual */
  residual(residual_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  /* print_scalar3d(residual_scalar3d); */
  gridtofourier(residual_coeff, residual_scalar3d, 0, 0);
  transform_fouriertoylm(residual_coeff, residual_ylm_coeff, fouriertoylm_matrix);
/*   printf("residual_ylm_coeff is here:\n"); */
/*   print_ylm_coeff(residual_ylm_coeff); */
  jacobian1(j1_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  jacobian2(j2_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  jacobian3(j3_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);


  /* evaluate a and effective source  */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  alpha = gsl_vector_get(alpha_vector, z);
    	  j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	  j2 = scalar3d_get(j2_scalar3d, z, i, j, k);
	  j3 = scalar3d_get(j3_scalar3d, z, i, j, k);
	  
	  a = alpha*alpha * (1.0 + j2*j2 + j3*j3) / (j1*j1);
	  
	  residual_d = scalar3d_get(residual_scalar3d, z, i, j, k);
	  source_d = scalar3d_get(source_scalar3d, z, i, j, k);
	  
	  roru = scalar3d_get(roru_scalar3d, z, i, j, k);
	  
	  if(z==nz-1&&i==nr-1) { /* infinity */
	    source_eff = (1.0*source_d + residual_d) / a;
	  } else if(z==nz-1) { /* external compactified zone */
	    xi_i = -cos(PI*i/(nr-1));
	    source_eff = (source_d/pow(roru, 4) + residual_d) * pow(alpha*(xi_i-1.0), 4) / a;
	  } else { /* kernel and shells */
	    source_eff = (source_d + residual_d) / a;
	  }
	  
	  scalar3d_set(source_eff_scalar3d, z, i, j, k, source_eff);
	  /* printf("z=%d, i=%d, j=%d, k=%d, roru=%.18e, alpha=%.18e, a=%.18e, res=%.18e, s=%.18e, s_eff=%.18e\n", z, i, j, k, roru, alpha, a, residual_d, source_d, source_eff); */
	  scalar3d_set(a_scalar3d, z, i, j, k, a);	
	}
      }
    }
  }

  gridtofourier(a_coeff, a_scalar3d, 0, 0);
  transform_fouriertoylm(a_coeff, a_ylm_coeff, fouriertoylm_matrix);
/*   printf("a_ylm_coeff is here:\n"); */
/*   print_ylm_coeff(a_ylm_coeff); */
  
  pseudo_spherical_laplacian(pseudolapf_scalar3d, field_scalar3d, alpha_vector, beta_vector);

  /* compare left and right hand sides */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  source_eff = scalar3d_get(source_eff_scalar3d, z, i, j, k);
 	  pseudolapf = scalar3d_get(pseudolapf_scalar3d, z, i, j, k);
	  error = (source_eff - pseudolapf)/source_eff;
	  /* printf("z=%d, i=%d, j=%d, k=%d, s_eff=%.18e, lapf_a=%.18e, err=%.18e\n", z, i, j, k, source_eff, pseudolapf, error); */
	}
      }
    }
  }
  
  /* solves \tilde \Delta f = \tilde \sigma where \tilde \sigma = (\sigma + R(f))/a */
  poisson_iteration(field_solution_scalar3d, source_eff_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* compare numerical to analytical solution */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  source_eff = scalar3d_get(source_scalar3d, z, i, j, k);
	  residual_d = scalar3d_get(residual_scalar3d, z, i, j, k);
 	  num = scalar3d_get(field_solution_scalar3d, z, i, j, k);
 	  anal = scalar3d_get(field_scalar3d, z, i, j, k);
	  error = (num - anal)/anal;
	  /* printf("z=%d, i=%d, j=%d, k=%d, xi_i=%.18e, theta_j=%.18e, phi_k=%.18e, s_eff=%.18e, res=%.18e, f_n=%.18e, f_a=%.18e, err=%.18e\n", z, i, j, k, xi_i, theta_j, phi_k, source_eff, residual_d, num, anal, error);*/
	}
      }
    }
  }

  /* compare numerical to analytical solution in terms of coefficients */
  gridtofourier(field_coeff, field_scalar3d, 0, 0);
  gridtofourier(field_solution_coeff, field_solution_scalar3d, 0, 0);
  k = 0;
  imag = 0;
  for (z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	num = coeff_get(field_solution_coeff, z, i, j, k, imag);
	anal = coeff_get(field_coeff, z, i, j, k, imag);
	error = (num - anal)/anal;
	/* printf("z=%d, i=%d, j=%d, k=%d, imag=%d, f_n=%.18e, f_a=%.18e, err=%.18e\n", z, i, j, k, imag, num, anal, error); */
      }
    }
  }

  /* compare numerical to analytical solution in terms of spherical harmonics */
  gridtofourier(field_coeff, field_scalar3d, 0, 0);
  gridtofourier(field_solution_coeff, field_solution_scalar3d, 0, 0);
  transform_fouriertoylm(field_coeff, field_ylm_coeff, fouriertoylm_matrix);
  transform_fouriertoylm(field_solution_coeff, field_solution_ylm_coeff, fouriertoylm_matrix);
  m = 0;
  imag = 0;
  /* for ( L = m; L < nt-m%2; L++ ) { */
  for ( L = m; L < 3; L++ ) {
    for ( z = 0; z < nz; z++ ) {
      for ( i = 0; i < nr; i++ ) {
	num = ylm_coeff_get(field_solution_ylm_coeff, z, i, L, m, imag);
	anal = ylm_coeff_get(field_ylm_coeff, z, i, L, m, imag);
	error = (num - anal)/anal;
	/* printf("z=%d, i=%d, L=%d, m=%d, imag=%d, f_n=%.18e, f_a=%.18e, err=%.18e\n", z, i, L, m, imag, num, anal, error); */
      }
    }
  }
  
  scalar2d_free(boundary_scalar2d);
  scalar2d_free(f_scalar2d);
  scalar2d_free(g_scalar2d);
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);
  scalar3d_free(field_scalar3d);
  scalar3d_free(source_scalar3d);
  scalar3d_free(roru_scalar3d);

  scalar3d_free(residual_scalar3d);
  scalar3d_free(j1_scalar3d);
  scalar3d_free(j2_scalar3d);
  scalar3d_free(j3_scalar3d);
  scalar3d_free(source_eff_scalar3d);
  scalar3d_free(pseudolapf_scalar3d);
  scalar3d_free(field_solution_scalar3d);
  coeff_free(field_solution_coeff);

  bound_coeff_free(f_bound_coeff);
  bound_coeff_free(g_bound_coeff);
  coeff_free(source_coeff);
  coeff_free(field_coeff);

  coeff_free(residual_coeff);
  ylm_coeff_free(residual_ylm_coeff);
  scalar3d_free(a_scalar3d);
  coeff_free(a_coeff);
  ylm_coeff_free(a_ylm_coeff);

  ylm_coeff_free(field_ylm_coeff);
  ylm_coeff_free(field_solution_ylm_coeff);

  fouriertoylm_matrix_free(fouriertoylm_matrix);
  ylmtofourier_matrix_free(ylmtofourier_matrix);

  return 0;
}

void pseudo_spherical_laplacian(scalar3d *pseudolapf_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector)
{
  int z, i, j, k;
  int nz, nr, nt, np;
  
  scalar3d *df_dxi_scalar3d;
  coeff *f_coeff;
  coeff *d2f_dxi2_coeff;
  scalar3d *d2f_dxi2_scalar3d;
  coeff *anglapf_coeff;
  scalar3d *anglapf_scalar3d;
  
  double alpha, beta, xi_i;

  double df_dxi;
  double d2f_dxi2;
  double anglapf;
  double lapf;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;

  df_dxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  f_coeff = coeff_alloc(nz, nr, nt, np);
  d2f_dxi2_coeff = coeff_alloc(nz, nr, nt, np);
  d2f_dxi2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  anglapf_coeff = coeff_alloc(nz, nr, nt, np);
  anglapf_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* calculate first derivative */
  dfdxi_ongrid(df_dxi_scalar3d, f_scalar3d);


  /* calculate second derivative */
  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  d2fdxi2(d2f_dxi2_coeff, f_coeff);
  fouriertogrid(d2f_dxi2_scalar3d, d2f_dxi2_coeff, 0, 0);

  /* calculate angular Laplacian */
  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  laplace_ang(anglapf_coeff, f_coeff);
  fouriertogrid(anglapf_scalar3d, anglapf_coeff, 0, 0);
  
  /* combine the terms */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if (z==0) {
	    alpha = gsl_vector_get(alpha_vector, z);
	    df_dxi = scalar3d_get(df_dxi_scalar3d, z, i, j, k); 
	    d2f_dxi2 = scalar3d_get(d2f_dxi2_scalar3d, z, i, j, k); 
	    anglapf = scalar3d_get(anglapf_scalar3d, z, i, j, k); 
	 
	    lapf = ( d2f_dxi2 + ( 2.0 / xi_i ) * df_dxi + (1.0 / pow(xi_i, 2)) * anglapf ) / (alpha*alpha);
	    scalar3d_set(pseudolapf_scalar3d, z, i, j, k, lapf);	  
	  } else if (z==nz-1) {
	    alpha = gsl_vector_get(alpha_vector, z);
	    d2f_dxi2 = scalar3d_get(d2f_dxi2_scalar3d, z, i, j, k); 
	    anglapf = scalar3d_get(anglapf_scalar3d, z, i, j, k); 
	    
	    /* lapf = ( d2f_dxi2 + (1.0 / pow((xi_i-1.0), 2)) * anglapf ) / (alpha*alpha); */
	    lapf = alpha*alpha * ( pow(xi_i-1.0, 4)*d2f_dxi2 + pow(xi_i-1.0, 2)*anglapf );
	    scalar3d_set(pseudolapf_scalar3d, z, i, j, k, lapf);	 
	  } else {
	    alpha = gsl_vector_get(alpha_vector, z);
	    beta = gsl_vector_get(beta_vector, z);
	    df_dxi = scalar3d_get(df_dxi_scalar3d, z, i, j, k); 
	    d2f_dxi2 = scalar3d_get(d2f_dxi2_scalar3d, z, i, j, k); 
	    anglapf = scalar3d_get(anglapf_scalar3d, z, i, j, k); 
	    
	    lapf = ( d2f_dxi2 + ( 2.0 / (xi_i+beta/alpha) ) * df_dxi + (1.0 / pow((xi_i+beta/alpha), 2)) * anglapf ) / (alpha*alpha);
	    scalar3d_set(pseudolapf_scalar3d, z, i, j, k, lapf);	  
	  }
	  
	}
      }
    }
  }
  
  scalar3d_free(df_dxi_scalar3d);
  coeff_free(f_coeff);
  coeff_free(d2f_dxi2_coeff);
  scalar3d_free(d2f_dxi2_scalar3d);
  coeff_free(anglapf_coeff);
  scalar3d_free(anglapf_scalar3d);
}

/**************************************************/
/* Take 2 partial derivatives with respect to xi. */
/*           d^2 f / d xi^2                       */
/**************************************************/
void d2fdxi2(coeff *d2fdxi2, coeff *f)
{
  int z;
  int i;
  int j;
  int k;
  int imag; /* imag=0 is real.  imag=1 is imaginary. */
  int nz;
  int nr; 
  int nt;
  int np;
  int npc;
  double fprime_i; /* f'_i(xi) */
  double f_ip1; /* f_{i+1}(xi) */
  double fprime_ip2; /* f'_{i+3}(xi) */
  
  coeff *dfdxi;

  nz = f->nz;
  nr = f->nr;
  nt = f->nt;
  np = f->np;
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  dfdxi = coeff_alloc(nz, nr, nt, np);

  /* even kernel */
  /* even Chebyshev polynomials become odd Chebyshev polynomials */
  /* even Chebyshev polynomials are indexed as T_{2i}(xi) */
  /* odd Chebyshev polynomials are indexed as T_{2i+1}(xi) */
  z = 0;
  for(imag=0; imag<=1; imag++) {
    for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */
      for(j = 0; j < nt; j+=2) { /* even j only */
	/* do recursion for xi derivatives */

	/* set dfdxi coefficient for T_{2(nr-1)+1} to zero: */
	coeff_set(dfdxi, z, nr-1, j, k, imag, 0.0);
	/* set dfdxi coefficient for T_{2(nr-2)}: */
	fprime_i = 4*(nr-1)*coeff_get(f, z, nr-1, j, k, imag);
	coeff_set(dfdxi, z, nr-2, j, k, imag, fprime_i);
	/* now set dfdxi for nr-3 and below: */
	fprime_ip2 = 0.0; /* where i starts at nr-3 below */
	for(i=nr-3; i>=0; i--) {
	  f_ip1 = coeff_get(f, z, i+1, j, k, imag);
	  fprime_ip2 = coeff_get(dfdxi, z, i+1, j, k, imag);
	  fprime_i = 4*(i+1)*f_ip1 + fprime_ip2;
	  coeff_set(dfdxi, z, i, j, k, imag, fprime_i);
	} 

	/* set dfdxi coefficient for T_{2(nr-1)} to zero: */
	coeff_set(d2fdxi2, z, nr-1, j, k, imag, 0.0);
	/* set dfdxi coefficient for T_{2(nr-2)+1}: */
	fprime_i = 2*(2*nr-3)*coeff_get(dfdxi, z, nr-2, j, k, imag);
	coeff_set(d2fdxi2, z, nr-2, j, k, imag, fprime_i);
	/* now set dfdxi for nr-3 and below: */
	fprime_ip2 = 0.0; /* where i starts at nr-3 below */
	for(i=nr-3; i>=0; i--) {
	  f_ip1 = coeff_get(dfdxi, z, i, j, k, imag);
	  fprime_ip2 = coeff_get(d2fdxi2, z, i+1, j, k, imag);
	  fprime_i = (2*(2*i+1)*f_ip1 + fprime_ip2)/(1+delta(i, 0));
	  coeff_set(d2fdxi2, z, i, j, k, imag, fprime_i);
	}

      }
    }
  }
  
  /* odd kernel */
  /* odd Chebyshev polynomials become even Chebyshev polynomials */
  for(imag=0; imag<=1; imag++) {
    for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */
      for(j = 1; j < nt-1; j+=2) { /* odd j only */
	/* do recursion for xi derivatives */
	/* set dfdxi coefficient for T_{2(nr-1)} to zero: */
	coeff_set(dfdxi, z, nr-1, j, k, imag, 0.0);
	/* set dfdxi coefficient for T_{2(nr-2)+1}: */
	fprime_i = 2*(2*nr-3)*coeff_get(f, z, nr-2, j, k, imag);
	coeff_set(dfdxi, z, nr-2, j, k, imag, fprime_i);
	/* now set dfdxi for nr-3 and below: */
	fprime_ip2 = 0.0; /* where i starts at nr-3 below */
	for(i=nr-3; i>=0; i--) {
	  f_ip1 = coeff_get(f, z, i, j, k, imag);
	  fprime_ip2 = coeff_get(dfdxi, z, i+1, j, k, imag);
	  fprime_i = (2*(2*i+1)*f_ip1 + fprime_ip2)/(1+delta(i, 0));
	  coeff_set(dfdxi, z, i, j, k, imag, fprime_i);
	} 

	/* set dfdxi coefficient for T_{2(nr-1)+1} to zero: */
	coeff_set(d2fdxi2, z, nr-1, j, k, imag, 0.0);
	/* set dfdxi coefficient for T_{2(nr-2)}: */
	fprime_i = 4*(nr-1)*coeff_get(dfdxi, z, nr-1, j, k, imag);
	coeff_set(d2fdxi2, z, nr-2, j, k, imag, fprime_i);
	/* now set dfdxi for nr-3 and below: */
	fprime_ip2 = 0.0; /* where i starts at nr-3 below */
	for(i=nr-3; i>=0; i--) {
	  f_ip1 = coeff_get(dfdxi, z, i+1, j, k, imag);
	  fprime_ip2 = coeff_get(d2fdxi2, z, i+1, j, k, imag);
	  fprime_i = 4*(i+1)*f_ip1 + fprime_ip2;
	  coeff_set(d2fdxi2, z, i, j, k, imag, fprime_i);
	} 

      }
    }
  }
  
  
  /* shells and external domain */
  for(z=1; z<nz; z++) {
    for(imag=0; imag<=1; imag++) {
      for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j = 0; j < nt; j++) { /* all j */
	  /* do recursion for xi derivatives */
	
	  /* set dfdxi for nr-1: */
	  coeff_set(dfdxi, z, nr-1, j, k, imag, 0.0);
	  /* set dfdxi for nr-2: */
	  fprime_i = 2*(nr-1)*coeff_get(f, z, nr-1, j, k, imag);
	  coeff_set(dfdxi, z, nr-2, j, k, imag, fprime_i);
	  /* now set dfdxi for nr-3 and below: */
	  fprime_ip2 = 0.0; /* where i starts at nr-3 below, so f'_{i+2} = f'_{nr-1} = 0 */
	  for(i=nr-3; i>=0; i--) {
	    f_ip1 = coeff_get(f, z, i+1, j, k, imag);
	    fprime_ip2 = coeff_get(dfdxi, z, i+2, j, k, imag);
	    fprime_i = (2*(i+1)*f_ip1 + fprime_ip2)/(1+delta(i, 0));
	    coeff_set(dfdxi, z, i, j, k, imag, fprime_i);
	  } 

	  /* set dfdxi for nr-1: */
	  coeff_set(d2fdxi2, z, nr-1, j, k, imag, 0.0);
	  /* set dfdxi for nr-2: */
	  fprime_i = 2*(nr-1)*coeff_get(dfdxi, z, nr-1, j, k, imag);
	  coeff_set(d2fdxi2, z, nr-2, j, k, imag, fprime_i);
	  /* now set dfdxi for nr-3 and below: */
	  fprime_ip2 = 0.0; /* where i starts at nr-3 below, so f'_{i+2} = f'_{nr-1} = 0 */
	  for(i=nr-3; i>=0; i--) {
	    f_ip1 = coeff_get(dfdxi, z, i+1, j, k, imag);
	    fprime_ip2 = coeff_get(d2fdxi2, z, i+2, j, k, imag);
	    fprime_i = (2*(i+1)*f_ip1 + fprime_ip2)/(1+delta(i, 0));
	    coeff_set(d2fdxi2, z, i, j, k, imag, fprime_i);
	  } 

	}
      }
    }
  }
  
  coeff_free(dfdxi);
}

/******************************/
/* Function for the boundary. */
/******************************/
/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 5.0; */
/*   else */
/*     return 10.0; */
/* } */
/* double boundary(int z, double theta, double phi) */
/* { */
/*   int L1 = 1; */
/*   int m1 = 0; */

/*   if(z==0) */
/*     return 5.0 + cos(theta); */
/*   else */
/*     return 10.0; */
/* } */
/* double boundary(int z, double theta, double phi) */
/* { */
/*   int L1 = 2; */
/*   int m1 = 1; */
  
/*   if(z==0) */
/*     return 5.0*(1.0 + 0.2*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))); */
/*   else if (z==1) */
/*     return 10.0; */
/*   else  */
/*     return 20.0*(1.0 + 0.2*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))); */
/* } */

double boundary(int z, double theta, double phi)
{
  if(z==0) {
    return 1.0;
  } else if (z==1) {
    return 5.0 + cos(theta);
  } else {
    return 10.0;
  }
}
/* double boundary(int z, double theta, double phi) */
/* { */
/*   int L1 = 1;  */
/*   int m1 = 0;  */
  
/*   if(z==0) { */
/*     return 1.0*(1.0 + 0.2*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))); */
/*   } else if (z==1) { */
/*     return 5.0; */
/*   } else if (z==2) { */
/*     return 10.0; */
/*   } else { */
/*     return 20.0; */
/*   } */
/* } */
/* double boundary(int z, double theta, double phi) */
/* { */
/*   double a; */
/*   double b; */
/*   double c; */
/*   double den; */

/*   phi = phi+1.0; */

/*   if(z==0) { */
/*     a = 1.7; */
/*     b = 1.0; */
/*     c = 1.2; */
/*     den = pow(sin(theta)*cos(phi) / a, 2) + pow(sin(theta)*sin(phi) / b, 2) + pow(cos(theta) / c, 2); */
/*     return pow(den, -0.5)*(1.0+0.5*cos(theta)); */
/*   } else if (z==1) { */
/*     a = 4.0; */
/*     b = 5.0; */
/*     c = 6.0; */
/*     den = pow(sin(theta)*cos(phi) / a, 2) + pow(sin(theta)*sin(phi) / b, 2) + pow(cos(theta) / c, 2); */
/*     return pow(den, -0.5)*(1.0+0.5*cos(theta)); */
/*   } else if (z==2) { */
/*     return 10.0; */
/*   } else { */
/*     a = 20.0; */
/*     b = 18.0; */
/*     c = 21.0; */
/*     den = pow(sin(theta)*cos(phi) / a, 2) + pow(sin(theta)*sin(phi) / b, 2) + pow(cos(theta) / c, 2); */
/*     return pow(den, -0.5)*(1.0+0.5*cos(theta)); */
/*   } */
/* } */

/**********************/
/* Source function.   */
/**********************/
double source(int z, double r, double theta, double phi)
{
  double R = 10.0;
  double rho = 1.0;
  if(z<3)
    return rho;
  else
    return 0.0;
}


/* double source(int z, double r, double theta, double phi) */
/* { */
/*   double R = 10.0; */
/*   if(z<3) */
/*     return (R - r*r/R); */
/*   else */
/*     return pow(R, 5)/pow(r, 4); */
/* } */
/* double source(int z, double r, double theta, double phi) */
/* { */
/*   int L1 = 3; */
/*   int m1 = 3; */
/*   int L2 = 4; */
/*   int m2 = 1; */
/*   double R = 10.0; */
/*   if(z<3) */
/*     return pow(r, L1)*((2*L1+3)*(2*L1+5)/pow(R, 2*L1+3) - (4*L1+10)*(2*L1+3)*r*r/pow(R, 2*L1+5)) */
/*       *gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)) */
/*       + pow(r, L2)*((2*L2+3)*(2*L2+5)/pow(R, 2*L2+3) - (4*L2+10)*(2*L2+3)*r*r/pow(R, 2*L2+5)) */
/*       *gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi)); */
/*   else */
/*     return 0.0; */
/* } */


/**********************************/
/* Solution to Poisson equation.  */
/**********************************/
double field(int z, double r, double theta, double phi)
{
  double R = 10.0;
  double rho = 1.0;
  if(z<3)
    return rho*(-R*R/2.0 + r*r/6.0);
  else
    return -rho*R*R*R/(3.0*r);
}

/* double field(int z, double r, double theta, double phi) */
/* { */
/*   double R = 10.0; */
/*   if(z<3) */
/*     return R*r*r/6.0 - pow(r, 4)/(20.0*R) - 3.0*pow(R, 3)/4.0; */
/*   else */
/*     return pow(R, 5)/(2.0*r*r) - 17.0*pow(R, 4)/(15.0*r); */
/* } */
/* double field(int z, double r, double theta, double phi) */
/* { */
/*   int L1 = 3; */
/*   int m1 = 3; */
/*   int L2 = 4; */
/*   int m2 = 1; */
/*   double R = 10.0; */
/*   if(z<3) */
/*     return pow(r, L1)*(0.5*(2*L1+5)*r*r/pow(R, 2*L1+3) - 0.5*(2*L1+3)*pow(r, 4)/pow(R, 2*L1+5)) */
/*       *gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)) */
/*       + pow(r, L2)*(0.5*(2*L2+5)*r*r/pow(R, 2*L2+3) - 0.5*(2*L2+3)*pow(r, 4)/pow(R, 2*L2+5)) */
/*       *gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi)); */
/*   else */
/*     return (1.0/pow(r, L1+1))*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)) */
/*       + (1.0/pow(r, L2+1))*gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi)); */
/* } */
