/**************************************************************
 * onebysintheta_dfdphi_test.c:                               *
 * Tests dividebysin, dfdphiprime, and fouriertogrid.         *
 *                                                            *
 * Author: Benjamin D. Lackey                                 *
 **************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h onebysintheta_dfdphi_test.c */

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
double onebysintheta_dfdphi(int z, double xi, double theta, double phi);

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
  coeff *onebysintheta_dfdphi_coeff;
  scalar3d *onebysintheta_dfdphi_scalar3d;
  
  double num, anal, error;
  
  nt = np/2 + 1;
  
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);
  dfdphi_coeff = coeff_alloc(nz, nr, nt, np);
  onebysintheta_dfdphi_coeff = coeff_alloc(nz, nr, nt, np);
  onebysintheta_dfdphi_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* evaluate field at gridpoints */
  functiontogrid_xi(field_scalar3d, field);

  /* determine coefficients of f */
  gridtofourier(field_coeff, field_scalar3d, 0, 0);
  printf("Original coefficients:\n");
  print_coeff(field_coeff);

  /* calculate df/dphi' */
  dfdphiprime(dfdphi_coeff, field_coeff);

  /* divide by sin(theta) */
  dividebysin(onebysintheta_dfdphi_coeff, dfdphi_coeff);
  printf("Coefficients after using 1/sin(theta) * df/dphi' operator:\n");
  print_coeff(onebysintheta_dfdphi_coeff);

  /* go back to gridpoints */
  fouriertogrid(onebysintheta_dfdphi_scalar3d, onebysintheta_dfdphi_coeff, 1, 1);
  
  /* compare numerical to analytic values of 1/sin(theta) * df/dphi' */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
 	  num = scalar3d_get(onebysintheta_dfdphi_scalar3d, z, i, j, k);
 	  anal = onebysintheta_dfdphi(z, xi_i, theta_j, phi_k);
	  error = (num - anal)/(anal);
	  printf("z=%d, i=%d, j=%d, k=%d, xi_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, xi_i, theta_j, phi_k, num, anal, error);
	}
      }
    }
  }
  
  scalar3d_free(field_scalar3d);
  coeff_free(field_coeff);
  coeff_free(dfdphi_coeff);
  coeff_free(onebysintheta_dfdphi_coeff);
  scalar3d_free(onebysintheta_dfdphi_scalar3d);
  
  return 0;
}

double field(int z, double xi, double theta, double phi)
{
  /*return xi*sin(theta)*(cos(phi) + sin(phi));*/
/*   return (cos(2.0*phi) + sin(2.0*phi)) */
/*     *(15.0/16.0)*(3.0 + 4.0*cos(2.0*theta) - 7.0*cos(4.0*theta)) /\* P^4_2(cos(theta)) *\/ */
/*     *xi*xi*xi*xi; */
  return (cos(3.0*phi) + sin(3.0*phi))
    *(-105.0/8.0)*(2.0*sin(2.0*theta) - sin(4.0*theta)) /* P^4_3(cos(theta)) */
    *xi;
}

double onebysintheta_dfdphi(int z, double xi, double theta, double phi)
{
  /*return xi*(-sin(phi) + cos(phi));*/
/*   return (-2.0*sin(2.0*phi) + 2.0*cos(2.0*phi)) */
/*     *(15.0/8.0)*(3.0*sin(theta) + 7.0*sin(3.0*theta)) /\* P^4_2(cos(theta)) / sin(theta) *\/ */
/*     *xi; */
  return (-3.0*sin(3.0*phi) + 3.0*cos(3.0*phi))
    *(-105.0/4.0)*(cos(theta) - cos(3.0*theta)) /* P^4_3(cos(theta)) / sin(theta) */
    *xi;
}
