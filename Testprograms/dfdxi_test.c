/**************************************************************************************
 * dfdxi_test.c:                                                                      *
 * Tests R/(xi+beta/alpha) by multiplying it by (xi+beta/alpha) and                   *
 * comparing it to R.                                                                 *
 *                                                                                    *
 * Author: Benjamin D. Lackey                                                         *
 **************************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/residual.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/poisson.h dfdxi_test.c */

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
double df_dxi(int z, double xi, double theta, double phi);
double onebyr_df_dxi_anal(int z, double xi, double theta, double phi, double alpha0);

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
  scalar3d *field_scalar3d;
  scalar3d *df_dxi_scalar3d;
  scalar3d *onebyr_df_dxi_scalar3d;
 
  coeff *r_coeff;
  coeff *dr_dxi_coeff;
  coeff *d2r_dxi2_coeff;
  scalar3d *d2r_dxi2_scalar3d;
  scalar3d *d2r_dxi2_anal_scalar3d;

  double alpha0;
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
  df_dxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_df_dxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  r_coeff = coeff_alloc(nz, nr, nt, np);
  dr_dxi_coeff = coeff_alloc(nz, nr, nt, np);
  d2r_dxi2_coeff = coeff_alloc(nz, nr, nt, np);
  d2r_dxi2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  d2r_dxi2_anal_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /*>>>>>>>>>>>> evaluate error of df_dx_ongrid function <<<<<<<<<<<<*/

  /* evaluate field(xi, theta, phi) at gridpoints */
  functiontogrid_xi(field_scalar3d, field);

/*   /\* evaluate df/dxi *\/ */
/*   dfdxi_ongrid(df_dxi_scalar3d, field_scalar3d); */

/*   /\* compare numerical to analytical value of df/dxi *\/ */
/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1))); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/*  	  num = scalar3d_get(df_dxi_scalar3d, z, i, j, k); */
/*  	  anal = df_dxi(z, xi_i, theta_j, phi_k); */
/* 	  error = (num - anal)/anal; */
/* 	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, xi_i, theta_j, phi_k, num, anal, error); */
/* 	} */
/*       } */
/*     } */
/*   } */


 /*  /\*>>>>>>>>>>>> evaluate error of onebyr_df_dxi function <<<<<<<<<<<<*\/ */

/*   /\* evaluate boundary function on gridpoints *\/ */
/*   boundarytogrid(boundary_scalar2d, boundary); */
  
/*   /\* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d *\/ */
/*   map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */
  
/*   /\* evaluate (1/R) * df/dxi *\/ */
/*   onebyr_df_dxi(onebyr_df_dxi_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */

/*   /\* now compare numerical to analytical value of (1/R) * df/dxi *\/ */
/*   rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */
/*   alpha0 = gsl_vector_get(alpha_vector, 0); */
/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  roru_i = scalar3d_get(r_scalar3d, z, i, j, k); */
/* 	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1))); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/*  	  num = scalar3d_get(onebyr_df_dxi_scalar3d, z, i, j, k); */
/*  	  anal = ((z==0 && i==0) ? onebyr_df_dxi_anal(z, xi_i, theta_j, phi_k, alpha0) : df_dxi(z, xi_i, theta_j, phi_k)/(roru_i)); */
/* 	  error = (num - anal)/anal; */
/* 	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, roru_i, theta_j, phi_k, num, anal, error); */
/* 	} */
/*       } */
/*     } */
/*   } */

  /*>>>>>>>>>>>> evaluate error of d2R_dxi2 function <<<<<<<<<<<<*/

  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* evaluate coefficients of R */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  gridtofourier(r_coeff, r_scalar3d, 0, 0);
  
  dfdxi(dr_dxi_coeff, r_coeff);
  dfdxi(d2r_dxi2_coeff, dr_dxi_coeff);
  
  /* go back to gridpoints */
  fouriertogrid(d2r_dxi2_scalar3d, d2r_dxi2_coeff, 0, 0);
  
  /* evaluate d^2R/dxi^2 with the analytical function */
  d2rdxi2(d2r_dxi2_anal_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
 	  num = scalar3d_get(d2r_dxi2_scalar3d, z, i, j, k);
 	  anal = scalar3d_get(d2r_dxi2_anal_scalar3d, z, i, j, k);;
	  error = (num - anal)/anal;
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
  scalar3d_free(df_dxi_scalar3d);
  scalar3d_free(onebyr_df_dxi_scalar3d);

  coeff_free(r_coeff);
  coeff_free(dr_dxi_coeff);
  coeff_free(d2r_dxi2_coeff);
  scalar3d_free(d2r_dxi2_scalar3d);
  scalar3d_free(d2r_dxi2_anal_scalar3d);

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

double field(int z, double xi, double theta, double phi)
{
  int l1 = 2;
  int m1 = 1;
  int l2 = 2;
  int m2 = 0;
  
  return pow(xi, l1)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
    + pow(xi, l2)*gsl_sf_legendre_sphPlm(l2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi));
}


double df_dxi(int z, double xi, double theta, double phi)
{
  int l1 = 2;
  int m1 = 1;
  int l2 = 2;
  int m2 = 0;

  return l1*pow(xi, l1-1)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
    + l2*pow(xi, l2-1)*gsl_sf_legendre_sphPlm(l2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi));
}

double onebyr_df_dxi_anal(int z, double xi, double theta, double phi, double alpha0)
{
  int l1 = 2;
  int m1 = 1;
  int l2 = 2;
  int m2 = 0;
  
  if(z==0 && xi<0.01)
    return l1*pow(xi, l1-2)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))/alpha0
      + l2*pow(xi, l2-2)*gsl_sf_legendre_sphPlm(l2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi))/alpha0;
  else  
    return l1*pow(xi, l1-1)*gsl_sf_legendre_sphPlm(l1, m1, cos(theta))*(cos(m1*phi) + sin(m1*phi))
      + l2*pow(xi, l2-1)*gsl_sf_legendre_sphPlm(l2, m2, cos(theta))*(cos(m2*phi) + sin(m2*phi));
  
  
 
}
