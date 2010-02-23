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


double boundary(int z, double theta, double phi);
double source(int z, double r, double theta, double phi);
double field(int z, double r, double theta, double phi);

int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 25; /* must be odd? */
  int nt;
  int np = 16; /* must be even */

  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *field_scalar3d;
  scalar3d *source_scalar3d;

  scalar3d *residual_scalar3d;
  scalar3d *j1_scalar3d;
  scalar3d *j2_scalar3d;
  scalar3d *j3_scalar3d;
  scalar3d *source_eff_scalar3d;
  scalar3d *field_solution_scalar3d;

  double alpha, j1, j2, j3, a, residual_d, source_d, source_eff;
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

  residual_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j3_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  source_eff_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_solution_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate source at gridpoints */
  functiontogrid(source_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, source);

  /* evaluate analytical solution at gridpoints */
  functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, field);

  /* evaluate residual */
  residual(residual_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
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
	  /*printf("z=%d, i=%d, j=%d, k=%d, alpha=%.18e\n", z, i, j, k, alpha);*/

	  residual_d = scalar3d_get(residual_scalar3d, z, i, j, k);
	  source_d = scalar3d_get(source_scalar3d, z, i, j, k);
	  
	  source_eff = (source_d + residual_d) / a;

	  scalar3d_set(source_eff_scalar3d, z, i, j, k, source_eff);
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
	  source_d = scalar3d_get(source_scalar3d, z, i, j, k);
 	  num = scalar3d_get(field_solution_scalar3d, z, i, j, k);
 	  anal = scalar3d_get(field_scalar3d, z, i, j, k);
	  error = (num - anal)/anal;
	  printf("z=%d, i=%d, j=%d, k=%d, s=%.18e, f_n=%.18e, f_a=%.18e, err=%.18e\n", z, i, j, k, source_d, num, anal, error);
	}
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

  scalar3d_free(residual_scalar3d);
  scalar3d_free(j1_scalar3d);
  scalar3d_free(j2_scalar3d);
  scalar3d_free(j3_scalar3d);
  scalar3d_free(source_eff_scalar3d);
  scalar3d_free(field_solution_scalar3d);

  return 0;
}



/******************************/
/* Function for the boundary. */
/******************************/
double boundary(int z, double theta, double phi)
{
  if(z==0)
    return 5.0;
  else
    return 10.0;
}
/* double boundary(int z, double theta, double phi) */
/* { */
/*   int L1 = 2; */
/*   int m1 = 1; */

/*   if(z==0) */
/*     return 1.0; */
/*   else if(z==1) */
/*     return 5.0/\**(1.0 + 0.7*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)))*\/; */
/*   else if(z==2) */
/*     return 10.0; */
/*   else */
/*     return 20.0; */
/* } */


/**********************/
/* Source function.   */
/**********************/
/* double source(int z, double r, double theta, double phi) */
/* { */
/*   double R = 10.0; */
/*   if(z<3) */
/*     return (R - r*r/R); */
/*   else */
/*     return pow(R, 5)/pow(r, 4); */
/* } */

double source(int z, double r, double theta, double phi)
{
  int L1 = 0;
  int m1 = 0;
  int L2 = 8;
  int m2 = 8;
  double R = 10.0;
  if(z<2)
    return pow(r, L1)*((2*L1+3)*(2*L1+5)/pow(R, 2*L1+3) - (4*L1+10)*(2*L1+3)*r*r/pow(R, 2*L1+5))
      *gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
      + pow(r, L2)*((2*L2+3)*(2*L2+5)/pow(R, 2*L2+3) - (4*L2+10)*(2*L2+3)*r*r/pow(R, 2*L2+5))
      *gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi));
  else
    return 0.0;
}


/**********************************/
/* Solution to Poisson equation.  */
/**********************************/
/* double field(int z, double r, double theta, double phi) */
/* { */
/*   double R = 10.0; */
/*   if(z<3) */
/*     return R*r*r/6.0 - pow(r, 4)/(20.0*R) - 3.0*pow(R, 3)/4.0; */
/*   else */
/*     return pow(R, 5)/(2.0*r*r) - 17.0*pow(R, 4)/(15.0*r); */
/* } */

double field(int z, double r, double theta, double phi)
{
  int L1 = 0;
  int m1 = 0;
  int L2 = 8;
  int m2 = 8;
  double R = 10.0;
  if(z<2)
    return pow(r, L1)*(0.5*(2*L1+5)*r*r/pow(R, 2*L1+3) - 0.5*(2*L1+3)*pow(r, 4)/pow(R, 2*L1+5))
      *gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
      + pow(r, L2)*(0.5*(2*L2+5)*r*r/pow(R, 2*L2+3) - 0.5*(2*L2+3)*pow(r, 4)/pow(R, 2*L2+5))
      *gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi));
  else
    return (1.0/pow(r, L1+1))*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
      + (1.0/pow(r, L2+1))*gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi));
}
