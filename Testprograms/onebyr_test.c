/**************************************************************
 * onebyr_test.c:                                             *
 * Tests function that divides by r after                     *
 * all types of basis function shifts.                        *
 *                                                            *
 * Author: Benjamin D. Lackey                                 *
 **************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h onebyr_test.c */

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

double field(int z, double r, double theta, double phi);
double fieldbyr(int z, double r, double theta, double phi);


int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 35; /* must be odd? */
  int nt;
  int np = 30; /* must be even */
  double r_i, theta_j, phi_k;
  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *field_scalar3d;
  coeff *field_coeff;
  scalar3d *fbyr_scalar3d;
  scalar3d *r_scalar3d;
  
  double num, anal, error;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);
  fbyr_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /*>>>>>>>>>>>>>>>>>> TESTING F/R ON NON-SHIFTED BASIS <<<<<<<<<<<<<<<<<<*/
  
  /* evaluate field at surface matched gridpoints */
  functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, field);

  /* find coefficients */
  gridtofourier(field_coeff, field_scalar3d, 0, 0);

  /* calculate f/R */
  dividebyr(fbyr_scalar3d, field_coeff, 0, 0, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  /* compare numerical to analytic values of the gradient at the surface matched gridpoints */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  num = scalar3d_get(fbyr_scalar3d, z, i, j, k);
	  anal = fieldbyr(z, r_i, theta_j, phi_k);
	  error = (num - anal)/anal;
	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, r_i, theta_j, phi_k, num, anal, error);
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
  coeff_free(field_coeff);
  scalar3d_free(fbyr_scalar3d);
  scalar3d_free(r_scalar3d);

  return 0;
}



/******************************/
/* Function for the boundary. */
/******************************/
/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0; */
/*   else */
/*     return 5.0; */
/* } */
double boundary(int z, double theta, double phi)
{
  if(z==0)
    return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/* + 0.2*cos(5*theta)*/);
  else
    return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/* + 0.2*cos(5*theta)*/);
}

/*>>>>>>>>>>> no basis shift <<<<<<<<<<<<<<*/

double field(int z, double r, double theta, double phi)
{
  if(z<2)
    return r*exp(-0.1*r*r)*cos(theta)*(cos(2*phi) + sin(2*phi));
  else
    return 1.0/(r*r);
}
double fieldbyr(int z, double r, double theta, double phi)
{
  if(z<2)
    return exp(-0.1*r*r)*cos(theta)*(cos(2*phi) + sin(2*phi));
  else
    return 1.0/r;
}
