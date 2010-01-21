/**************************************************************
 * gradient.c:                                                *
 * Functions to calculate the gradient of a scalar in         *
 * spherical coordinates.                                     *
 *                                                            *
 * Author: Benjamin D. Lackey                                 *
 **************************************************************/

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


/**************************************************/
/* Evaluates r component of the gradient of f     */
/* at the surface matched gridpoints.             */
/**************************************************/
void gradient_r(scalar3d *gradf_r_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k;
  int nz, nr, nt, np;
  coeff *f_coeff;
  coeff *dfdxi_coeff;
  scalar3d *dfdxi_scalar3d;
  scalar3d *j1_scalar3d;
  scalar3d *roru_scalar3d;
  double u, dfdxi_double, j1;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;

  f_coeff = coeff_alloc(nz, nr, nt, np);
  dfdxi_coeff = coeff_alloc(nz, nr, nt, np);
  dfdxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  roru_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  dfdxi(dfdxi_coeff, f_coeff);
  fouriertogrid(dfdxi_scalar3d, dfdxi_coeff, 1, 0);
  jacobian1(j1_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  /*print_scalar3d(j1_scalar3d);*/
  
  /* set df/dr = J_1^{-1} * df/dxi except for external domain */
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  dfdxi_double = scalar3d_get(dfdxi_scalar3d, z, i, j, k);
	  j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	  scalar3d_set(gradf_r_scalar3d, z, i, j, k, dfdxi_double/j1);
	}
      }
    }
  }
 
  /* set df/dr = - U^2 * J_1^{-1} * df/dxi in the external domain */
  rofxtp(roru_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  z = nz-1;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	u = scalar3d_get(roru_scalar3d, z, i, j, k);
	dfdxi_double = scalar3d_get(dfdxi_scalar3d, z, i, j, k);
	j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	scalar3d_set(gradf_r_scalar3d, z, i, j, k, -u*u*dfdxi_double/j1);
      }
    }
  }

  coeff_free(f_coeff);
  coeff_free(dfdxi_coeff);
  scalar3d_free(dfdxi_scalar3d);
  scalar3d_free(j1_scalar3d);
  scalar3d_free(roru_scalar3d);
}


/**************************************************/
/* Evaluates theta component of the gradient of f */
/* at the surface matched gridpoints.             */
/**************************************************/
void gradient_theta(scalar3d *gradf_theta_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k;
  int nz, nr, nt, np;
  coeff *f_coeff;
  coeff *dfdxi_coeff;
  scalar3d *dfdxi_scalar3d;
  coeff *dfdt_coeff;
  scalar3d *dfdt_byr_scalar3d;
  scalar3d *j1_scalar3d;
  scalar3d *j2_scalar3d;
  scalar3d *roru_scalar3d;
  double u, dfdt_byr, dfdxi_double, j1, j2;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;

  f_coeff = coeff_alloc(nz, nr, nt, np);
  dfdxi_coeff = coeff_alloc(nz, nr, nt, np);
  dfdxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  dfdt_coeff = coeff_alloc(nz, nr, nt, np);
  dfdt_byr_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  roru_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* calculate df/dxi by going to coefficients then back to gridpoints */
  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  dfdxi(dfdxi_coeff, f_coeff);
  fouriertogrid(dfdxi_scalar3d, dfdxi_coeff, 1, 0);
  
  /* calculate (1/R)*df/dtheta' */
  dfdthetaprime(dfdt_coeff, f_coeff);
  dividebyr(dfdt_byr_scalar3d, dfdt_coeff, 0, 1, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* calculate the xi and theta part of the Jacobian */
  jacobian1(j1_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  /*print_scalar3d(j1_scalar3d);*/
  jacobian2(j2_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  /*print_scalar3d(j2_scalar3d);*/

  /* set df/dtheta = (1/R)*df/dtheta' - (J_2/J_1)*df/dxi except for external domain */
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  dfdt_byr = scalar3d_get(dfdt_byr_scalar3d, z, i, j, k);
	  dfdxi_double = scalar3d_get(dfdxi_scalar3d, z, i, j, k);
	  j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	  j2 = scalar3d_get(j2_scalar3d, z, i, j, k);
	  scalar3d_set(gradf_theta_scalar3d, z, i, j, k, dfdt_byr - dfdxi_double*j2/j1);
	}
      }
    }
  }
 
  /* set df/dtheta = U^2 * ((1/R)*df/dtheta' - (J_2/J_1)*df/dxi) in the external domain */
  rofxtp(roru_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  z = nz-1;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	u = scalar3d_get(roru_scalar3d, z, i, j, k);
	dfdt_byr = scalar3d_get(dfdt_byr_scalar3d, z, i, j, k);
	dfdxi_double = scalar3d_get(dfdxi_scalar3d, z, i, j, k);
	j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	j2 = scalar3d_get(j2_scalar3d, z, i, j, k);
	scalar3d_set(gradf_theta_scalar3d, z, i, j, k, u*u*(dfdt_byr - dfdxi_double*j2/j1));
      }
    }
  }

  coeff_free(f_coeff);
  coeff_free(dfdxi_coeff);
  scalar3d_free(dfdxi_scalar3d);
  coeff_free(dfdt_coeff);
  scalar3d_free(dfdt_byr_scalar3d);
  scalar3d_free(j1_scalar3d);
  scalar3d_free(j2_scalar3d);
  scalar3d_free(roru_scalar3d);
}


/**************************************************/
/* Evaluates phi component of the gradient of f   */
/* at the surface matched gridpoints.             */
/**************************************************/
void gradient_phi(scalar3d *gradf_phi_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k;
  int nz, nr, nt, np;
  coeff *f_coeff;
  coeff *dfdxi_coeff;
  scalar3d *dfdxi_scalar3d;
  coeff *dfdp_coeff;
  coeff *dfdp_bysint_coeff;
  scalar3d *dfdp_byrsint_scalar3d;
  scalar3d *j1_scalar3d;
  scalar3d *j3_scalar3d;
  scalar3d *roru_scalar3d;
  double u, dfdp_byrsint, dfdxi_double, j1, j3;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;

  f_coeff = coeff_alloc(nz, nr, nt, np);
  dfdxi_coeff = coeff_alloc(nz, nr, nt, np);
  dfdxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  dfdp_coeff = coeff_alloc(nz, nr, nt, np);
  dfdp_bysint_coeff = coeff_alloc(nz, nr, nt, np);
  dfdp_byrsint_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j3_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  roru_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* calculate df/dxi by going to coefficients then back to gridpoints */
  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  /*print_coeff(f_coeff);*/
  dfdxi(dfdxi_coeff, f_coeff);
  fouriertogrid(dfdxi_scalar3d, dfdxi_coeff, 1, 0);
  
  /* calculate 1/(R*sin(theta)) * df/dphi' */
  dfdphiprime(dfdp_coeff, f_coeff);
  dividebysin(dfdp_bysint_coeff, dfdp_coeff);
  /*print_coeff(dfdp_bysint_coeff);*/
  dividebyr(dfdp_byrsint_scalar3d, dfdp_bysint_coeff, 1, 1, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* calculate the xi and theta part of the Jacobian */
  jacobian1(j1_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  jacobian3(j3_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* set df/dtheta = 1/(R*sin(theta)) * df/dphi' - (J_2/J_1)*df/dxi except for external domain */
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  dfdp_byrsint = scalar3d_get(dfdp_byrsint_scalar3d, z, i, j, k);
	  dfdxi_double = scalar3d_get(dfdxi_scalar3d, z, i, j, k);
	  j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	  j3 = scalar3d_get(j3_scalar3d, z, i, j, k);
	  scalar3d_set(gradf_phi_scalar3d, z, i, j, k, dfdp_byrsint - dfdxi_double*j3/j1);
	}
      }
    }
  }
 
  /* set df/dtheta = U^2 * 1/(R*sin(theta)) * df/dphi' - (J_2/J_1)*df/dxi) in the external domain */
  rofxtp(roru_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  z = nz-1;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	u = scalar3d_get(roru_scalar3d, z, i, j, k);
	dfdp_byrsint = scalar3d_get(dfdp_byrsint_scalar3d, z, i, j, k);
	dfdxi_double = scalar3d_get(dfdxi_scalar3d, z, i, j, k);
	j1 = scalar3d_get(j1_scalar3d, z, i, j, k);
	j3 = scalar3d_get(j3_scalar3d, z, i, j, k);
	scalar3d_set(gradf_phi_scalar3d, z, i, j, k, u*u*(dfdp_byrsint - dfdxi_double*j3/j1));
      }
    }
  }

  coeff_free(f_coeff);
  coeff_free(dfdxi_coeff);
  scalar3d_free(dfdxi_scalar3d);
  coeff_free(dfdp_coeff);
  coeff_free(dfdp_bysint_coeff);
  scalar3d_free(dfdp_byrsint_scalar3d);
  scalar3d_free(j1_scalar3d);
  scalar3d_free(j3_scalar3d);
  scalar3d_free(roru_scalar3d);
}
