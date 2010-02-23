/**************************************************************************************
 * r_xiplusbbya_test.c:                                                               *
 * Tests R/(xi+beta/alpha) by multiplying it by (xi+beta/alpha) and                   *
 * comparing it to R.                                                                 *
 *                                                                                    *
 * Author: Benjamin D. Lackey                                                         *
 **************************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/residual.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/poisson.h r_xiplusbbya_test.c */

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

int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 25; /* must be odd? */
  int nt;
  int np = 16; /* must be even */
  double roru_i, xi_i, theta_j, phi_k;
  
  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *r_scalar3d;
  scalar3d *r_xiplusbbya_scalar3d;
  
  double alpha, beta;
  double rfromfunction, error;
  
  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  r_xiplusbbya_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate R at gridpoints */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate R/(xi+beta/alpha) */
  r_xiplusb_a(r_xiplusbbya_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* compare functions with f=R */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  roru_i = scalar3d_get(r_scalar3d, z, i, j, k);
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  if(z==0) {
	    rfromfunction = scalar3d_get(r_xiplusbbya_scalar3d, z, i, j, k)*(xi_i);
	  } else if(z==nz-1) {
	    rfromfunction = scalar3d_get(r_xiplusbbya_scalar3d, z, i, j, k)*(xi_i - 1.0);
	  } else {
	    alpha = gsl_vector_get(alpha_vector, z);
	    beta = gsl_vector_get(beta_vector, z);
	    rfromfunction = scalar3d_get(r_xiplusbbya_scalar3d, z, i, j, k)*(xi_i + beta/alpha);
	  }
	  
	  error = (rfromfunction - roru_i)/roru_i;
	  printf("z=%d, i=%d, j=%d, k=%d, rfunc=%.18e, r=%.18e, err=%.18e\n", z, i, j, k, rfromfunction, roru_i, error);
	}
      }
    }
  }
  
  
  scalar2d_free(boundary_scalar2d);
  scalar2d_free(f_scalar2d);
  scalar2d_free(g_scalar2d);
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);
  scalar3d_free(r_scalar3d);
  scalar3d_free(r_xiplusbbya_scalar3d);
  
  return 0;
}


double boundary(int z, double theta, double phi)
{
  int L1 = 2;
  int m1 = 1;
  int L2 = 5;
  int m2 = 0;
  
  if(z==0)
    return 5.0*(1.0 
		+ 0.2*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
		+ 0.2*gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi)));
  else
    return 10.0*(1.0 
		+ 0.2*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
		+ 0.2*gsl_sf_legendre_sphPlm(L2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi)));
}
