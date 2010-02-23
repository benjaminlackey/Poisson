/*********************************************************************************
 * Makes sure that the mapping reproduces the boundary.                          *
 *                                                                               *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/residual.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/poisson.h mapsurface_test.c */

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
  int z, j, k;
  int nz = 3;
  int nt;
  int np = 16; /* must be even */

  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  
  double f, g, alpha, beta;
  double in, out, anal, error;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  
  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  print_vector(alpha_vector); 
  print_vector(beta_vector);
  print_scalar2d(f_scalar2d); 
  print_scalar2d(g_scalar2d);

  /* compare numerical to analytical solution */ 
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      /* bound of kernel */
      z = 0;
      f = scalar2d_get(f_scalar2d, z, j, k);
      g = scalar2d_get(g_scalar2d, z, j, k);
      alpha = gsl_vector_get(alpha_vector, z);
      
      out = alpha*(1.0 + f + g);
      anal = scalar2d_get(boundary_scalar2d, z, j, k);
      error = (out - anal)/anal;
      printf("z=%d, j=%d, k=%d, out_anal=%.18e, out_num=%.18e, err=%.18e\n", z, j, k, out, anal, error);
      
      /* bounds of shells */
      for ( z = 1; z < nz-1; z++ ) {
	  f = scalar2d_get(f_scalar2d, z, j, k);
	  g = scalar2d_get(g_scalar2d, z, j, k);
	  alpha = gsl_vector_get(alpha_vector, z);
	  beta = gsl_vector_get(beta_vector, z);

	  in = alpha*(-1.0 + f) + beta;
	  anal = scalar2d_get(boundary_scalar2d, z-1, j, k);
	  error = (in - anal)/anal;
	  printf("z=%d, j=%d, k=%d, in_anal=%.18e, in_num=%.18e, err=%.18e\n", z, j, k, in, anal, error);

	  out = alpha*(1.0 + g) + beta;
	  anal = scalar2d_get(boundary_scalar2d, z, j, k);
	  error = (out - anal)/anal;
	  printf("z=%d, j=%d, k=%d, out_anal=%.18e, out_num=%.18e, err=%.18e\n", z, j, k, out, anal, error);
      }

      /* inner bound of external domain */
      z = nz-1;
      f = scalar2d_get(f_scalar2d, z, j, k);
      alpha = gsl_vector_get(alpha_vector, z);
	
      in = 1.0 / (alpha*(-2.0 + f));
      anal = scalar2d_get(boundary_scalar2d, z-1, j, k);
      error = (in - anal)/anal;
      printf("z=%d, j=%d, k=%d, in_anal=%.18e, in_num=%.18e, err=%.18e\n", z, j, k, in, anal, error);
    }
  }
  

  scalar2d_free(boundary_scalar2d);
  scalar2d_free(f_scalar2d);
  scalar2d_free(g_scalar2d);
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);

  return 0;
}


double boundary(int z, double theta, double phi)
{
  int L1 = 3;
  int m1 = 2;
  if(z==0)
    return 5.0*(1.0 + 0.3*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)));
  else
    return 10.0*(1.0 + 0.3*gsl_sf_legendre_sphPlm(L1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)));
}
