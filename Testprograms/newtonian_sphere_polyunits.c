/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/poisson.h newtoniansphericaltest.c */

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

double boundary(int z, int nt, int np, double theta, double phi);
double rho_func(int z, double r, double theta, double phi);
void functiontogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi));

int main (void)
{
  
  FILE *fpgrid;
  
  int z, i, j, k;
  int nz = 2;
  int nr = 25; /* must be odd? */
  int nt;
  int np = 12; /* must be even */
  double r_i, theta_j, phi_k;
  scalar2d *surface_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *r_scalar3d;
  scalar3d *rho_scalar3d;
  scalar3d *source_scalar3d;
  scalar3d *phi_scalar3d;
  scalar3d *enthalpy_scalar3d;
  coeff *source_coeff;
  coeff *phi_coeff;
  ylm_coeff *source_ylm_coeff;
  ylm_coeff *phi_ylm_coeff;
  gsl_matrix **fouriertoylm;
  gsl_matrix **ylmtofourier;

  int iteration;
  double c;
  double rho_max;
  double enthalpy_max;
  double enthalpy;
  double rho_new;
  double polyindex = 1;

  int rsurface;
  int zsurface;
  double enthalpy_absmin;
  
  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  surface_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  rho_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  source_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  phi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  enthalpy_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  source_coeff = coeff_alloc(nz, nr, nt, np);
  phi_coeff = coeff_alloc(nz, nr, nt, np);
  source_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  phi_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);

  /* make matrices for converting between double fourier series and spherical harmonics */ 
  fouriertoylm = fouriertoylm_matrix_alloc(nt);
  ylmtofourier = ylmtofourier_matrix_alloc(nt);
  fouriertoylm_matrix_set(fouriertoylm);
  ylmtofourier_matrix_set(ylmtofourier);
  
  /* Evaluate analytic boundary function on the boundary grid points. */
  for ( z = 0; z < nz-1; z++ ) {
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	scalar2d_set(surface_scalar2d, z, j, k, boundary(nt, np, z, theta_j, phi_k));
      }
    }
  }
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(surface_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  /*print_vector(alpha_vector);
  print_vector(beta_vector);
  print_scalar2d(f_scalar2d);
  print_scalar2d(g_scalar2d);*/

  /***********************************************************/
  /* give initial guess for density and radius */
  /* fix rho_max forever */
  /* start of iteration */
  /*   calculate phi from density */
  /*   phi at surface gives C (H=0 there) using current guess of radius from density profile */
  /*   phi and C then give H everywhere */
  /*   new guess for surface defined by H=0 */
  /*   H gives rho from EOS */
  /* end of iteration */
  /***********************************************************************/

  /*>>>>>>>>>>>>>>>>> FIX CENTRAL DENSITY AND SET INITIAL GUESSES <<<<<<<<<<<<<<<<<<*/

  /* fill grid with initial guess for density profile */
  functiontogrid(rho_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, rho_func);
  /*print_scalar3d(source_scalar3d);*/

  /* calculate initial source guess S = 4*PI*rho */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(source_scalar3d, z, i, j, k, 
		       4*PI*scalar3d_get(rho_scalar3d, z, i, j, k));
	  }
      }
    }
  }  
  
  rho_max = 1.0;
  enthalpy_max = pow(rho_max, 1/polyindex);

  zsurface=1;
  rsurface=0;

  /*>>>>>>>>>>>>>>> BEGIN ITERATION <<<<<<<<<<<<<<<<<<*/

  for ( iteration = 0; iteration < 25; iteration++ ) {
    
    /* rewrite function in terms of coefficients of Chebyshev polynomials and double Fourier series */
    gridtofourier(source_coeff, source_scalar3d, 0, 0);
    /*print_coeff(source_coeff);*/
 
    /* go from fourier series to spherical harmonics */
    transform_fouriertoylm(source_coeff, source_ylm_coeff, fouriertoylm);
    /*print_ylm_coeff(source_ylm_coeff);*/
    
    /* solve poisson equation and return coefficients */
    solve_poisson_spherical(phi_ylm_coeff, source_ylm_coeff, alpha_vector, beta_vector);
    
    /* go from spherical harmonics to fourier series */
    transform_ylmtofourier(phi_ylm_coeff, phi_coeff, ylmtofourier);
    
    /* fourier and Chebyshev coefficients --> values on grid */
    fouriertogrid(phi_scalar3d, phi_coeff, 0, 0);

    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    /*HERE IS THE PROBLEM.  THE INNER SURFACE OF THE 2ND DOMAIN IS NOT NECESSARILY THE SURFACE OF THE STAR.*/
    /* just outside the surface at arbitrary {theta, phi} */ 
    c = scalar3d_get(phi_scalar3d, zsurface, rsurface, 2, 2); 
      
    /* calculate enthalpy H = C - Phi */
    for ( z = 0; z < nz; z++ ) {
      for ( i = 0; i < nr; i++ ) {
	for ( j = 0; j < nt; j++ ) {
	  for ( k = 0; k < np; k++ ) {
	    scalar3d_set(enthalpy_scalar3d, z, i, j, k, 
			 c - scalar3d_get(phi_scalar3d, z, i, j, k));
	  }
	}
      }
    }  
    /*print_scalar3d(enthalpy_scalar3d);*/
    
    /* determine zsurface, rsurface by the position where enthalpy is closest to zero */
    rsurface = 0;
    zsurface = 0;
    enthalpy_absmin = scalar3d_get(enthalpy_scalar3d, 0, 0, 2, 2);
    printf("\n");
    for ( z = 0; z < nz; z++ ) {
      for ( i = 0; i < nr; i++ ) {
	enthalpy = scalar3d_get(enthalpy_scalar3d, z, i, 2, 2);
	printf("%f, ", enthalpy);
	if(fabs(enthalpy) < fabs(enthalpy_absmin)) {
	  enthalpy_absmin = enthalpy;
	  zsurface = z;
	  rsurface = i;
	}
      }
    }
    printf("\n");
    printf("enthalpy=%f, zsurface=%d, rsurface=%d\n", enthalpy_absmin, zsurface, rsurface); 
    
    /* set new source term */
    for ( z = 0; z < nz; z++ ) {
      for ( i = 0; i < nr; i++ ) {
	for ( j = 0; j < nt; j++ ) {
	  for ( k = 0; k < np; k++ ) {
	    enthalpy = scalar3d_get(enthalpy_scalar3d, z, i, j, k);
	    if(enthalpy >= 0.0) { /* matter exists only where enthalpy is positive */
	      rho_new = pow(enthalpy/enthalpy_max, polyindex);
	      scalar3d_set(source_scalar3d, z, i, j, k, 4*PI*rho_new);
	    } else { /* no matter */
	      scalar3d_set(source_scalar3d, z, i, j, k, 0.0); /* source = 4*PI*0.0 */
	    }
	  }
	}
      }
    }  
    /*print_scalar3d(source_scalar3d);*/
    
  } 

/*>>>>>>>>>>>>>>>>>>>>>>>> end of iteration <<<<<<<<<<<<<<<<<<<<<<<<<<<*/ 
  
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* print to file the position and function value at all gridpoints */
  fpgrid=fopen("rhoofrtp.txt", "w");
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  r_i = scalar3d_get(r_scalar3d, z, i, j, k);
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  fprintf(fpgrid, "%.18e\t%.18e\t%.18e\t%.18e\n", r_i, theta_j, phi_k, scalar3d_get(source_scalar3d, z, i, j, k)/(4*PI));
	}
      }
    }
  }
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	r_i = 1.0 / scalar3d_get(r_scalar3d, z, i, j, k);
	theta_j = PI*j/(nt-1);
	phi_k = 2*PI*k/np;
	fprintf(fpgrid, "%.18e\t%.18e\t%.18e\t%.18e\n", r_i, theta_j, phi_k, scalar3d_get(source_scalar3d, z, i, j, k)/(4*PI));
      }
    }
  }
  
  /* free memory */ 
  /*   scalar2d_free(surface_scalar2d); */
  /*   scalar2d_free(f_scalar2d); */
  /*   scalar2d_free(g_scalar2d); */
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);
  scalar3d_free(r_scalar3d);
  scalar3d_free(rho_scalar3d);
  scalar3d_free(source_scalar3d);
  scalar3d_free(phi_scalar3d);
  scalar3d_free(enthalpy_scalar3d);
  coeff_free(source_coeff);
  coeff_free(phi_coeff);
  ylm_coeff_free(source_ylm_coeff);
  ylm_coeff_free(phi_ylm_coeff);
  
  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
double boundary(int nt, int np, int z, double theta, double phi)
{
  return 1.0;
}

/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(0==z) */
/*     return 5.0; */
/*   else */
/*     return 10.0; */
/* } */

/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(0==z) */
/*     return 3.0; */
/*   else if(1==z) */
/*     return 5.0; */
/*   else */
/*     return 10.0; */
/* } */


/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else if(z==1) */
/*     return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/* } */



/*********************************************************/
/* Some function of position in spherical coordinates.   */
/*********************************************************/
/* double rho_func(int z, double r, double theta, double phi) */
/* { */
/*   int L = 5; */
/*   double R = 5.0; */
/*   if(0==z) */
/*     return gsl_sf_legendre_sphPlm(L, 0, cos(theta))*0.5*(L+3)*(L+5)*((L-4)*r*r/pow(R, L+5) - (L-2)/pow(R, L+3)); */
/*   else */
/*     return 0.0; */
/* } */

/* double rho_func(int z, double r, double theta, double phi) */
/* { */
/*   if(0==z) */
/*     return sin(PI*(r+0.0001))/(PI*(r+0.0001)); */
/*   else */
/*     return 0.0; */
/* } */

double rho_func(int z, double r, double theta, double phi)
{
  if(0==z)
    return 1.0;
  else
    return 0.0;
}


/******************************************************************************************************/
/* Take a function f = f(z, r, theta, phi) and evaluate it on the surface matched grid func_scalar3d. */
/******************************************************************************************************/ 
void functiontogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi))
{
  int z, i, j, k;
  int nz, nr, nt, np;
  double xi_i, theta_j, phi_k;
  double r;
  scalar3d *r_scalar3d;

  nz = func_scalar3d->nz;
  nr = func_scalar3d->nr;
  nt = func_scalar3d->nt;
  np = func_scalar3d->np;
  
  /* value of R (or U) at each point (xi, theta, phi) on grid */
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);  
  
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* kernel and shells */
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if(z==0)
	    xi_i = sin(PI*i/(2*(nr-1)));
	  else
	    xi_i = -cos(PI*i/(nr-1));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  r = scalar3d_get(r_scalar3d, z, i, j, k);
	  scalar3d_set(func_scalar3d, z, i, j, k, func(z, r, theta_j, phi_k));
	}
      }
    }
  }
  /* external domain except for r = infinity */
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	xi_i = -cos(PI*i/(nr-1));
	theta_j = PI*j/(nt-1);
	phi_k = 2*PI*k/np;
	r = scalar3d_get(r_scalar3d, z, i, j, k);
	scalar3d_set(func_scalar3d, z, i, j, k, func(z, 1/r, theta_j, phi_k));
      }
    }
  }
  /* set to zero at r = infinity */
  z=nz-1;
  i=nr-1;
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      scalar3d_set(func_scalar3d, z, i, j, k, 0.0);
    }
  }

  scalar3d_free(r_scalar3d);
}
