/**************************************************************
 * dfdphi_test.c:                                             *
 * Tests dfdphiprime, and fouriertogrid.                      *
 *                                                            *
 * Author: Benjamin D. Lackey                                 *
 **************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h dfdphi_test.c */

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

double field(int z, double xi, double theta, double phi);
double dfield_dphi(int z, double xi, double theta, double phi);

int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 5; /* must be odd? */
  int nt;
  int np = 10; /* must be even */
  double xi_i, theta_j, phi_k;
  scalar3d *field_scalar3d;
  coeff *field_coeff;
  coeff *dfdphi_coeff;
  scalar3d *dfdphi_scalar3d;
  
  double dfdphi, dfdphi_analytic, error;
  
  nt = np/2 + 1;
  
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);
  dfdphi_coeff = coeff_alloc(nz, nr, nt, np);
  dfdphi_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* evaluate field at gridpoints */
  functiontogrid_xi(field_scalar3d, field);

  /* determine coefficients of f */
  gridtofourier(field_coeff, field_scalar3d, 0, 0);
  print_coeff(field_coeff);

  /* calculate df/dphi' */
  dfdphiprime(dfdphi_coeff, field_coeff);
  print_coeff(dfdphi_coeff);

  /* go back to gridpoints */
  fouriertogrid(dfdphi_scalar3d, dfdphi_coeff, 0, 0);
	  
 /* compare numerical to analytic values of df/dphi' */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
 	  dfdphi = scalar3d_get(dfdphi_scalar3d, z, i, j, k);
 	  dfdphi_analytic = dfield_dphi(z, xi_i, theta_j, phi_k);
	  error = (dfdphi - dfdphi_analytic)/(dfdphi_analytic);
	  printf("z=%d, i=%d, j=%d, k=%d, xi_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, xi_i, theta_j, phi_k, dfdphi, dfdphi_analytic, error);
	}
      }
    }
  }
  
  scalar3d_free(field_scalar3d);
  coeff_free(field_coeff);
  coeff_free(dfdphi_coeff);
  scalar3d_free(dfdphi_scalar3d);
  
  return 0;
}

double field(int z, double xi, double theta, double phi)
{
  /*return 1.0;*/
  return cos(2.0*phi);
  /*return cos(phi)*/
  /*return cos(theta);*/
  /*return xi*sin(theta)*(cos(phi) + sin(phi));*/
}

double dfield_dphi(int z, double xi, double theta, double phi)
{
  /*return 0.0;*/
  return -2.0*sin(2.0*phi);
  /*return -sin(phi);*/
  /*return 0.0;*/
  /*return xi*sin(theta)*(-sin(phi) + cos(phi));*/
}
