/***************************************************************************
 * Test the onebyrsq_anglaplacef function                                  *
 *   Should return 0 for r^0 Y_00                                          *
 *   Not defined for l=1                                                   *
 *   Should return well defined value for r^2 Y_2m                         *
 *   Should return 0 for r^l Y_lm with l>2                                 *
 *                                                                         *
 * Author: Benjamin D. Lackey                                              *
 ***************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h /Users/lackey/Research/Poisson/residual.c ang_lap_test.c */

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
double field(int z, double xi, double theta, double phi);
double onebyr2_lap_field(int z, double xi, double theta, double phi, double alpha0);

int main (void)
{
int z, i, j, k;
  int nz = 3;
  int nr = 25; /* must be odd? */
  int nt;
  int np = 16; /* must be even */
  double r_i, xi_i, theta_j, phi_k;
  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *r_scalar3d;
  scalar3d *field_scalar3d;
  scalar3d *onebyr2_lapfield_scalar3d;
  double alpha0, num, anal, error;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr2_lapfield_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate field(xi, theta, phi) at gridpoints */
  functiontogrid_xi(field_scalar3d, field);

  /* evaluate (1/R^2)*nabla_{theta' phi'} */
  onebyrsq_anglaplacef(onebyr2_lapfield_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* now compare numerical to analytical solution of (1/R^2)*nabla_{theta' phi'} */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);  
  alpha0 = gsl_vector_get(alpha_vector, 0);
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  r_i = scalar3d_get(r_scalar3d, z, i, j, k);
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
 	  num = scalar3d_get(onebyr2_lapfield_scalar3d, z, i, j, k);
 	  anal = ((z==0 && i==0) ? onebyr2_lap_field(z, xi_i, theta_j, phi_k, alpha0) : onebyr2_lap_field(z, xi_i, theta_j, phi_k, alpha0)/(r_i*r_i));
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
  scalar3d_free(r_scalar3d);
  scalar3d_free(field_scalar3d);
  scalar3d_free(onebyr2_lapfield_scalar3d);

  return 0;
}

/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 5.0; */
/*   else */
/*     return 10.0; */
/* } */

double boundary(int z, double theta, double phi)
{
  int l = 2;
  int m = 1;

  if(z==0)
    return 5.0*(1.0 + 0.1*gsl_sf_legendre_sphPlm(l, m, cos(theta))*(cos(m*phi) + sin(m*phi)));
  else
    return 10.0*(1.0 + 0.1*gsl_sf_legendre_sphPlm(l, m, cos(theta))*(cos(m*phi) + sin(m*phi)));
}

/* double field(int z, double xi, double theta, double phi) */
/* { */
/*   int l1 = 2; */
/*   int m1 = 0; */
/*   int l2 = 5; */
/*   int m2 = 4; */
  
/*   if(z<2) */
/*     return pow(r, l1)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)) */
/*       + pow(r, l2)*gsl_sf_legendre_sphPlm(l2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi)); */
/*   else */
/*     return (1/pow(r, l1))*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)) */
/*       + (1/pow(r, l2))*gsl_sf_legendre_sphPlm(l2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi)); */
/* } */

double field(int z, double xi, double theta, double phi)
{
  int l1 = 2;
  int m1 = 1;
  
  return pow(xi, l1)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi));
}


double onebyr2_lap_field(int z, double xi, double theta, double phi, double alpha0)
{
  int l1 = 2;
  int m1 = 1;
  
  if(z==0 && xi<0.01)
    return -l1*(l1+1.0)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi)) / (alpha0*alpha0);
  else  
    return -l1*(l1+1.0)*pow(xi, l1)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi));
}
