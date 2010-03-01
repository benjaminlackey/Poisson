/*********************************************************************************
 * Test the expression R(f) by evaluating both sides of the equation             *
 *             a \tilde\Delta f = \sigma(f) + R(f)                               *
 *                                                                               *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h /Users/lackey/Research/Poisson/residual.c residual_test.c */

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
  scalar3d *field_saved_scalar3d;
  scalar3d *residual_scalar3d;

  double num, anal, error;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_saved_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  residual_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate field and source at gridpoints */
  functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, field);
  functiontogrid(source_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, field);

  /* calculate remainder from field */
  residual(residual_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
 
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
 	  num = scalar3d_get(field_scalar3d, z, i, j, k);
 	  anal = scalar3d_get(field_saved_scalar3d, z, i, j, k);
	  error = (num - anal)/(anal);
	  printf("z=%d, i=%d, j=%d, k=%d, %.18e, %.18e, %.18e\n", z, i, j, k, num, anal, error);
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
  scalar3d_free(field_saved_scalar3d);
  scalar3d_free(residual_scalar3d);

  return 0;
}


/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0; */
/*   else */
/*     return 3.0; */
/* } */
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


double field(int z, double r, double theta, double phi)
{
  int l1 = 2;
  int m1 = 1;
  int l2 = 2;
  int m2 = 0;
  
  return pow(r, l1)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
    + pow(r, l2)*gsl_sf_legendre_sphPlm(l2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi));
}


laplacian_tilde(scalar3d *laptilde_scalar3d, scalar3d *field_scalar3d)
{
  
}
