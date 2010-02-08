/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/residual.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/poisson.h poissonsphericalboundtest.c */

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
double source_func(int z, double r, double theta, double phi);
double field(int z, double r, double theta, double phi);


int main (void)
{
  
  FILE *fpgrid;
  
  int z, i, j, k;
  int nz = 5;
  int nr = 15; /* must be odd? */
  int nt;
  int np = 16; /* must be even */
  double r_i, theta_j, phi_k;
  scalar2d *surface_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *r_scalar3d;
  scalar3d *source_scalar3d;
  scalar3d *field_scalar3d;
  coeff *field_coeff;
  coeff *source_coeff;
  ylm_coeff *field_ylm_coeff;
  ylm_coeff *source_ylm_coeff;
  gsl_matrix **fouriertoylm;
  gsl_matrix **ylmtofourier;

  double s, num, anal, error;
  
  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  surface_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  source_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  source_coeff = coeff_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);
  source_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);
  field_ylm_coeff = ylm_coeff_alloc(nz, nr, nt, np);

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
/*   print_vector(alpha_vector); */
/*   print_vector(beta_vector); */
/*   print_scalar2d(f_scalar2d); */
/*   print_scalar2d(g_scalar2d); */

  /* fill grid with data points determined by the function field */
  functiontogrid(source_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, source_func);
  /*print_scalar3d(source_scalar3d);*/

  /* rewrite function in terms of coefficients of Chebyshev polynomials and double Fourier series */
  gridtofourier(source_coeff, source_scalar3d, 0, 0);
  /* print_coeff(source_coeff); */

  /* go from fourier series to spherical harmonics */
  transform_fouriertoylm(source_coeff, source_ylm_coeff, fouriertoylm);
  /* print_ylm_coeff(source_ylm_coeff); */

  /* solve poisson equation and return coefficients */
  solve_poisson_spherical(field_ylm_coeff, source_ylm_coeff, alpha_vector, beta_vector);
  
  /* go from sphericas harmonics to fourier series */
  transform_ylmtofourier(field_ylm_coeff, field_coeff, ylmtofourier);
  
  /* fourier and Chebyshev coefficients --> values on grid */
  fouriertogrid(field_scalar3d, field_coeff, 0, 0);
  
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* compare numerical to analytical solution */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  s = scalar3d_get(field_scalar3d, z, i, j, k);
 	  num = scalar3d_get(field_scalar3d, z, i, j, k);
 	  anal = field(z, r_i, theta_j, phi_k);
	  error = (num - anal)/anal;
	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, r_i, theta_j, phi_k, num, anal, error);
	}
      }
    }
  }

/*   /\* print to file the position and function value at all gridpoints *\/ */
/*   fpgrid=fopen("fofrtp.txt", "w"); */
/*   for ( z = 0; z < nz-1; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  r_i = scalar3d_get(r_scalar3d, z, i, j, k); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/* 	  fprintf(fpgrid, "%.18e\t%.18e\t%.18e\t%.18e\n", r_i, theta_j, phi_k, scalar3d_get(field_scalar3d, z, i, j, k)); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   z=nz-1; */
/*   for ( i = 0; i < nr-1; i++ ) { */
/*     for ( j = 0; j < nt; j++ ) { */
/*       for ( k = 0; k < np; k++ ) { */
/* 	r_i = 1.0 / scalar3d_get(r_scalar3d, z, i, j, k); */
/* 	theta_j = PI*j/(nt-1); */
/* 	phi_k = 2*PI*k/np; */
/* 	fprintf(fpgrid, "%.18e\t%.18e\t%.18e\t%.18e\n", r_i, theta_j, phi_k, scalar3d_get(field_scalar3d, z, i, j, k)); */
/*       } */
/*     } */
/*   } */
  
  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*     return 5.0; */
/* } */

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
/*   else */
/*     return 5.0; */
/* } */

double boundary(int nt, int np, int z, double theta, double phi)
{
  if(0==z)
    return 1.0;
  else if(1==z)
    return 3.0;
  else if(2==z)
    return 5.0;
  else
    return 8.0;
}


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
/* double source_func(int z, double r, double theta, double phi) */
/* { */
/*   int L = 9; */
/*   double R = 5.0; */
/*   if(z<3) */
/*   else */
/*     return 0.0; */
/* } */

/* double source_func(int z, double r, double theta, double phi) */
/* { */
/*   int L = 4; */
/*   int m = 3; */

/*   return pow(r, L)*gsl_sf_legendre_sphPlm(L, m, cos(theta))*cos(m*phi); */
/* } */

double source_func(int z, double r, double theta, double phi)
{
  int L = 4;
  int m = 3;
  double R = 5.0;
  if(z<3)
    return pow(r, L)*((2*L+3)*(2*L+5)/pow(R, 2*L+3) - (4*L+10)*(2*L+3)*r*r/pow(R, 2*L+5))
      *gsl_sf_legendre_sphPlm(L, m, cos(theta))*cos(m*phi);
  else
    return 0.0;
}

/* double source_func(int z, double r, double theta, double phi) */
/* { */
/*   double R = 5.0; */
/*   if(0==z) */
/*     return (R - r*r/R); */
/*   else if(1==z) */
/*     return (R - r*r/R); */
/*   else */
/*     return pow(R, 5)/pow(r, 4); */
/* } */

/* double field(int z, double r, double theta, double phi) */
/* { */
/*   int L = 9; */
/*   double R = 5.0; */
/*   if(z<3) */
/*     return gsl_sf_legendre_sphPlm(L, 0, cos(theta))*((L+5)*r*r/(2.0*pow(R, L+3)) - (L+3)*pow(r, 4)/(2.0*pow(R, L+5))); */
/*   else */
/*     return gsl_sf_legendre_sphPlm(L, 0, cos(theta))/pow(r, L+1); */
/* } */

double field(int z, double r, double theta, double phi)
{
  int L = 4;
  int m = 3;
  double R = 5.0;
  if(z<3)
    return pow(r, L)*(0.5*(2*L+5)*r*r/pow(R, 2*L+3) - 0.5*(2*L+3)*pow(r, 4)/pow(R, 2*L+5))
      *gsl_sf_legendre_sphPlm(L, m, cos(theta))*cos(m*phi);
  else
    return (1.0/pow(r, L+1))*gsl_sf_legendre_sphPlm(L, m, cos(theta))*cos(m*phi);
}
