/*********************************************************************************
 * Test the laprace_ang_bound function                                           *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h /Users/lackey/Research/Poisson/residual.c anglaplace_bound_test.c */

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

double boundary(int z, double theta, double phi);
double lap_boundary(int z, double theta, double phi);

int main (void)
{
  int z, j, k;
  int nz = 3;
  int nt;
  int np = 16; /* must be even */
  double theta_j, phi_k;
  scalar2d *f_scalar2d;
  bound_coeff *f_bound_coeff;
  bound_coeff *lapf_bound_coeff;
  scalar2d *lapf_scalar2d;
  double num, anal, error;

  nt = np/2 + 1;
  
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  f_bound_coeff = bound_coeff_alloc(nz, nt, np);
  lapf_bound_coeff = bound_coeff_alloc(nz, nt, np);
  lapf_scalar2d = scalar2d_alloc(nz, nt, np);

  /* evaluate boundary function on gridpoints */
  boundarytogrid(f_scalar2d, boundary);
  
  /* convert to coefficients */
  gridtofourier_bound(f_bound_coeff, f_scalar2d);

  /* take angular Laplacian */
  laplace_ang_bound(lapf_bound_coeff, f_bound_coeff);

  /* go back to gridpoints */
  fouriertogrid_bound(lapf_scalar2d, lapf_bound_coeff, 0);

  /* compare numerical to analytic values */
  for ( z = 0; z < nz; z++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	theta_j = PI*j/(nt-1);
	phi_k = 2*PI*k/np;
	num = scalar2d_get(lapf_scalar2d, z, j, k);
	anal = lap_boundary(z, theta_j, phi_k);
	error = (num - anal)/(anal);
	printf("z=%d, j=%d, k=%d, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, j, k, theta_j, phi_k, num, anal, error);
      }
    }
  }

  scalar2d_free(f_scalar2d);
  bound_coeff_free(f_bound_coeff);
  bound_coeff_free(lapf_bound_coeff);
  scalar2d_free(lapf_scalar2d);

  return 0;
}


double boundary(int z, double theta, double phi)
{
  double x;
  int l = 5;
  int m = 3;
  
  x = cos(theta);
  return gsl_sf_legendre_sphPlm(l, m, x)*(cos(m*phi) + sin(m*phi));
}

double lap_boundary(int z, double theta, double phi)
{
  double x;
  int l = 5;
  int m = 3;
  
  x = cos(theta);
  return -l*(l+1.0)*gsl_sf_legendre_sphPlm(l, m, x)*(cos(m*phi) + sin(m*phi));
}
