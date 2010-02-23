/*********************************************************************************
 * onebyrsintheta_d2f_dphidxi_test.c:                                            *
 * Compares onebyrsin_d2rbydpdx to onebyrsintheta_d2f_dphidxi.                   *
 * Both functions return same values when R(xi, theta, phi) = f(xi, theta, phi). * 
 *                                                                               *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/residual.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/poisson.h onebyr_d2f_dthetadxi_test.c */

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
double d2f_dthetadxi(int z, double xi, double theta, double phi);

int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 35; /* must be odd? */
  int nt;
  int np = 30; /* must be even */
  double roru_i, xi_i, theta_j, phi_k;
  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *r_scalar3d;
  scalar3d *field_scalar3d;
  scalar3d *onebyr_d2f_dthetadxi_scalar3d;
  scalar3d *outf_scalar3d;
  scalar3d *outr_scalar3d;

  double num, anal, error;
  double outf, outr, diff;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_d2f_dthetadxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  outf_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  outr_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* evaluate field */
  functiontogrid_xi(field_scalar3d, field);

  /* evaluate the function R(xi, theta, phi */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

/*   /\*>>>>>>>>>>>> compare numerical to analytical <<<<<<<<<<<<<*\/ */
/*   onebyr_d2f_dthetadxi(onebyr_d2f_dthetadxi_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */

/*   /\* compare them *\/ */
/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  roru_i = scalar3d_get(r_scalar3d, z, i, j, k); */
/* 	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1))); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/* 	  num = scalar3d_get(onebyr_d2f_dthetadxi_scalar3d, z, i, j, k) * roru_i; */
/* 	  anal = d2f_dthetadxi(z, xi_i, theta_j, phi_k); */
/* 	  error = (num - anal)/anal; */
/* 	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n",  */
/* 		 z, i, j, k, roru_i, theta_j, phi_k, num, anal, error); */
/* 	} */
/*       } */
/*     } */
/*   } */

  /*>>>>> compare both versions of 1/R d^2R/(dtheta'dxi) with 2 different functions <<<<<<*/

  onebyr_d2f_dthetadxi(outf_scalar3d, r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyr_d2r_dthetadxi(outr_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* compare them */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  roru_i = scalar3d_get(r_scalar3d, z, i, j, k);
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  outf = scalar3d_get(outf_scalar3d, z, i, j, k);
	  outr = scalar3d_get(outr_scalar3d, z, i, j, k);
	  diff = (outf - outr)/outr;
	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, roru_i, theta_j, phi_k, outf, outr, diff);
	}
      }
    }
  }

  scalar2d_free(boundary_scalar2d);
  scalar2d_free(f_scalar2d);
  scalar2d_free(g_scalar2d);
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);
  scalar3d_free(onebyr_d2f_dthetadxi_scalar3d);
  scalar3d_free(r_scalar3d);
  scalar3d_free(field_scalar3d);
  scalar3d_free(outr_scalar3d);
  scalar3d_free(outf_scalar3d);

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
/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/\* + 0.2*cos(5*theta)*\/); */
/*   else */
/*     return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/\* + 0.2*cos(5*theta)*\/); */
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
/* double boundary(int z, double theta, double phi) */
/* { */
/*   double a; */
/*   double b; */
/*   double c; */
/*   double den; */

/*   if(z==0) { */
/*     a = 5.0; */
/*     b = 6.0; */
/*     c = 7.0; */
/*     den = pow(sin(theta)*cos(phi) / a, 2) + pow(sin(theta)*sin(phi) / b, 2) + pow(cos(theta) / c, 2); */
/*     return pow(den, -0.5); */
/*   } else { */
/*     a = 10.0; */
/*     b = 8.0; */
/*     c = 9.0; */
/*     den = pow(sin(theta)*cos(phi) / a, 2) + pow(sin(theta)*sin(phi) / b, 2) + pow(cos(theta) / c, 2); */
/*     return pow(den, -0.5); */
/*   } */
/* } */

/***************************************/
/* Field and operator acting on field. */
/***************************************/
double field(int z, double xi, double theta, double phi)
{
  int l1 = 2;
  int m1 = 2;
  int l2 = 1;
  int m2 = 1;
  
  return pow(xi, l1)*sin(theta)*sin(theta)*(cos(m1*phi) + sin(m1*phi)) /* prop x^2 Y_2^2 */
    + pow(xi, l2)*sin(theta)*(cos(m2*phi) + sin(m2*phi)); /* prop x^1 Y_1^0 */
}


double d2f_dthetadxi(int z, double xi, double theta, double phi)
{
  int l1 = 2;
  int m1 = 2;
  int l2 = 1;
  int m2 = 1;

  return l1*pow(xi, l1-1)*2.0*sin(theta)*cos(theta)*(cos(m1*phi) + sin(m1*phi))
    + l2*pow(xi, l2-1)*cos(theta)*(cos(m2*phi) + sin(m2*phi));
}
