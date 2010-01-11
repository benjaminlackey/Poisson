/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c  /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/poisson.h remapsurfacetest.c */

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

double boundary(int nt, int np, int z, double theta, double phi);
double enthalpy_func(int z, double r, double theta, double phi);
void functiontogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi));

int main (void)
{
  
  FILE *fpgrid;
  
  int z, j, k;
  int nz = 3;
  int nr = 75; /* must be odd? */
  int nt;
  int np = 12; /* must be even */
  double r_i, theta_j, phi_k;
  scalar2d *surface_scalar2d; 
  scalar2d *newsurface_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;

  scalar3d *enthalpy_scalar3d;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  surface_scalar2d = scalar2d_alloc(nz-1, nt, np);
  newsurface_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  enthalpy_scalar3d = scalar3d_alloc(nz, nr, nt, np);

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
  
  /* fill grid with initial guess for density profile */
  functiontogrid(enthalpy_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, enthalpy_func);
  
  /* find new surface */
  findnewsurface(newsurface_scalar2d, enthalpy_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* print to file the position and function value at all gridpoints */
  fpgrid=fopen("surfaceoftp.txt", "w");
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      z=0;
      theta_j = PI*j/(nt-1);
      phi_k = 2*PI*k/np;
      r_i = scalar2d_get(newsurface_scalar2d, z, j, k);
      fprintf(fpgrid, "%.18e\t%.18e\t%.18e\n", theta_j, phi_k, r_i);
    }
  }
  
  /* free memory */ 
  /*   scalar2d_free(surface_scalar2d); */
  /*   scalar2d_free(f_scalar2d); */
  /*   scalar2d_free(g_scalar2d); */
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);
  scalar3d_free(enthalpy_scalar3d);
  
  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/

double boundary(int nt, int np, int z, double theta, double phi)
{
  if(0==z)
    return 1.0e6;
  else
    return 1.0e7;
}


/*********************************************************/
/* Some function of position in spherical coordinates.   */
/*********************************************************/
double enthalpy_func(int z, double r, double theta, double phi)
{
  return 2.8e14*(1.2e6*(1.0 - 0.8*cos(theta)*cos(theta)) - r);
}


/******************************************************************************************************/
/* Take a function f = f(z, r, theta, phi) and evaluate it on the surface matched grid func_scalar3d. */
/******************************************************************************************************/ 
void functiontogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi))
{
  int z, i, j, k;
  int nz, nr, nt, np;
  double xi_i, theta_j, phi_k;
  double r;
  scalar3d *r_scalar3d;

  nz = func_scalar3d->nz;
  nr = func_scalar3d->nr;
  nt = func_scalar3d->nt;
  np = func_scalar3d->np;
  
  /* value of R (or U) at each point (xi, theta, phi) on grid */
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);  
  
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* kernel and shells */
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if(z==0)
	    xi_i = sin(PI*i/(2*(nr-1)));
	  else
	    xi_i = -cos(PI*i/(nr-1));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  r = scalar3d_get(r_scalar3d, z, i, j, k);
	  scalar3d_set(func_scalar3d, z, i, j, k, func(z, r, theta_j, phi_k));
	}
      }
    }
  }
  /* external domain except for r = infinity */
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	xi_i = -cos(PI*i/(nr-1));
	theta_j = PI*j/(nt-1);
	phi_k = 2*PI*k/np;
	r = scalar3d_get(r_scalar3d, z, i, j, k);
	scalar3d_set(func_scalar3d, z, i, j, k, func(z, 1/r, theta_j, phi_k));
      }
    }
  }
  /* set to zero at r = infinity */
  z=nz-1;
  i=nr-1;
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      scalar3d_set(func_scalar3d, z, i, j, k, 0.0);
    }
  }

  scalar3d_free(r_scalar3d);
}
