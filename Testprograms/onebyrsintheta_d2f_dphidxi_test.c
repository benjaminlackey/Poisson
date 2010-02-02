/*********************************************************************************
 * onebyrsintheta_d2f_dphidxi_test.c:                                            *
 * Compares onebyrsintheta_d2f_dphidxi to analytical function.                   *
 *                                                                               *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h onebyrsintheta_d2f_dphidxi_test.c */

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
double field(int z, double xi, double theta, double phi);
double onebysintheta_d2field_dphidxi(int z, double xi, double theta, double phi);


int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 35; /* must be odd? */
  int nt;
  int np = 30; /* must be even */
  double xi_i, theta_j, phi_k;
  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *r_scalar3d;
  scalar3d *field_scalar3d;
  scalar3d *out_scalar3d;
  scalar3d *onebysintheta_d2f_dphidxi_scalar3d;
  double num, anal, error;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  out_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebysintheta_d2f_dphidxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate field at gridpoints */
  functiontogrid_xi(field_scalar3d, field);

  /* take derivatives */
  onebyrsintheta_d2f_dphidxi(out_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate the function R(xi, theta, phi */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* multiply by R to get d^2f / (dtheta dxi) */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(onebysintheta_d2f_dphidxi_scalar3d, z, i, j, k,
		       scalar3d_get(r_scalar3d, z, i, j, k)*scalar3d_get(out_scalar3d, z, i, j, k));
	}
      }
    }
  }


  /* compare numerical to analytic values of d^2f/(dtheta' dxi) */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
 	  num = scalar3d_get(onebysintheta_d2f_dphidxi_scalar3d, z, i, j, k);
 	  anal = onebysintheta_d2field_dphidxi(z, xi_i, theta_j, phi_k);
	  error = (num - anal)/(anal);
	  printf("z=%d, i=%d, j=%d, k=%d, xi_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, xi_i, theta_j, phi_k, num, anal, error);
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
  scalar3d_free(field_scalar3d);
  scalar3d_free(out_scalar3d);
  scalar3d_free(onebysintheta_d2f_dphidxi_scalar3d);
  
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
/*     return 3.0; */
/* } */
double boundary(int z, double theta, double phi)
{
  if(z==0)
    return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/* + 0.2*cos(5*theta)*/);
  else
    return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/* + 0.2*cos(5*theta)*/);
}



/*********************************************/
/* Some function in spherical coordinates.   */
/*********************************************/
double field(int z, double xi, double theta, double phi)
{
  return xi*xi*sin(theta)*sin(theta)*(cos(2.0*phi) + sin(2.0*phi));
}

/**************************************/
/* The corresponding theta gradient.  */
/**************************************/
double onebysintheta_d2field_dphidxi(int z, double xi, double theta, double phi)
{
  return 2.0*xi*sin(theta)*(-2.0*sin(2.0*phi) + 2.0*cos(2.0*phi));
}
