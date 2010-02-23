/*********************************************************************************
 * Test the residual function by comparing it to Mathematica notebook            *
 * which calculates the residual analytically.                                   *
 *                                                                               *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h /Users/lackey/Research/Poisson/residual.c residual_test_type2.c */

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

double ffunc(int z, double theta, double phi);
double gfunc(int z, double theta, double phi);
double field(int z, double xi, double theta, double phi);

int main (void)
{
  int z, i, j, k; 
  int nz = 3;
  int nr = 55; /* must be odd? */
  int nt;
  int np = 30; /* must be even */

  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *field_scalar3d;

  scalar3d *j1_scalar3d;
  scalar3d *j2_scalar3d;
  scalar3d *j3_scalar3d;
  scalar3d *rby_xiplusbbya_scalar3d;
  scalar3d *df_dxi_scalar3d;
  scalar3d *onebyr_df_dxi_scalar3d;
  scalar3d *d2r_dxi2_scalar3d;
  scalar3d *onebyr_d2f_dthetadxi_scalar3d;
  scalar3d *onebyr_d2r_dthetadxi_scalar3d;
  scalar3d *onebyrsintheta_d2f_dphidxi_scalar3d;
  scalar3d *onebyrsintheta_d2r_dphidxi_scalar3d;
  scalar3d *onebyrsq_anglaplacef_scalar3d;
  scalar3d *onebyrsq_anglaplacer_scalar3d;

  scalar3d *residual_scalar3d;

  double roru_i, xi_i, theta_j, phi_k;
  double num, anal, error;

  nt = np/2 + 1;
  
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j3_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  rby_xiplusbbya_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  df_dxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_df_dxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  d2r_dxi2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_d2f_dthetadxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_d2r_dthetadxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsintheta_d2f_dphidxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsintheta_d2r_dphidxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsq_anglaplacef_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsq_anglaplacer_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  residual_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* evaluate boundary function on gridpoints */
  boundarytogrid(f_scalar2d, ffunc);
  boundarytogrid(g_scalar2d, gfunc);
  print_scalar2d(f_scalar2d);

  gsl_vector_set(alpha_vector, 0, 1.0);
  gsl_vector_set(alpha_vector, 1, 1.0);
  gsl_vector_set(alpha_vector, 2, -0.02);
  gsl_vector_set(beta_vector, 0, 0.0);
  gsl_vector_set(beta_vector, 1, 8.0);
  gsl_vector_set(beta_vector, 2, 0.0);
  
  /* evaluate field at gridpoints */
  functiontogrid_xi(field_scalar3d, field);

  /* calculate remainder from field */
  jacobian1(j1_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  jacobian2(j2_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  jacobian3(j3_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  r_xiplusb_a(rby_xiplusbbya_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  dfdxi_ongrid(df_dxi_scalar3d, field_scalar3d);
  onebyr_df_dxi(onebyr_df_dxi_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  d2rdxi2(d2r_dxi2_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  onebyr_d2f_dthetadxi(onebyr_d2f_dthetadxi_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyr_d2r_dthetadxi(onebyr_d2r_dthetadxi_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyrsintheta_d2f_dphidxi(onebyrsintheta_d2f_dphidxi_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyrsintheta_d2r_dphidxi(onebyrsintheta_d2r_dphidxi_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyrsq_anglaplacef(onebyrsq_anglaplacef_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyrsq_anglaplacer(onebyrsq_anglaplacer_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  residual(residual_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
 
  /*>>>>>>>>>> kernel <<<<<<<<<<<*/
  z = 0;
  i = 9;
  j = 5;
  k = 9;
  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
  theta_j = PI*j/(nt-1);
  phi_k = 2*PI*k/np;
  printf("z=%d, i=%d, j=%d, k=%d, x=%.18e, t=%.18e, p=%.18e\n", z, i, j, k, xi_i, theta_j, phi_k);
  printf("%.18e\n", scalar3d_get(j1_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(j2_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(j3_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(rby_xiplusbbya_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(df_dxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_df_dxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(d2r_dxi2_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_d2f_dthetadxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_d2r_dthetadxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsintheta_d2f_dphidxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsintheta_d2r_dphidxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsq_anglaplacef_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsq_anglaplacer_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(residual_scalar3d, z, i, j, k));

  /*>>>>>>>>>> shell <<<<<<<<<<<<*/
  z = 1;
  i = 9;
  j = 5;
  k = 9;
  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
  theta_j = PI*j/(nt-1);
  phi_k = 2*PI*k/np;
  printf("z=%d, i=%d, j=%d, k=%d, x=%.18e, t=%.18e, p=%.18e\n", z, i, j, k, xi_i, theta_j, phi_k);
  printf("%.18e\n", scalar3d_get(j1_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(j2_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(j3_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(rby_xiplusbbya_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(df_dxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_df_dxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(d2r_dxi2_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_d2f_dthetadxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_d2r_dthetadxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsintheta_d2f_dphidxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsintheta_d2r_dphidxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsq_anglaplacef_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsq_anglaplacer_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(residual_scalar3d, z, i, j, k));
  
  /*>>>>>>>>> external <<<<<<<<<<*/
  z = 2;
  i = 9;
  j = 5;
  k = 9;
  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
  theta_j = PI*j/(nt-1);
  phi_k = 2*PI*k/np;
  printf("z=%d, i=%d, j=%d, k=%d, x=%.18e, t=%.18e, p=%.18e\n", z, i, j, k, xi_i, theta_j, phi_k);
  printf("%.18e\n", scalar3d_get(j1_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(j2_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(j3_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(rby_xiplusbbya_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(df_dxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_df_dxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(d2r_dxi2_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_d2f_dthetadxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyr_d2r_dthetadxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsintheta_d2f_dphidxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsintheta_d2r_dphidxi_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsq_anglaplacef_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(onebyrsq_anglaplacer_scalar3d, z, i, j, k));
  printf("%.18e\n", scalar3d_get(residual_scalar3d, z, i, j, k));


/*   /\* compare numerical to analytic values *\/ */
/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1))); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/*  	  num = scalar3d_get(residual_scalar3d, z, i, j, k); */
/* 	  printf("z=%d, i=%d, j=%d, k=%d, x=%.18e, t=%.18e, p=%.18e, %.18e\n", z, i, j, k, xi_i, theta_j, phi_k, num); */
/* 	} */
/*       } */
/*     } */
/*   } */


  scalar2d_free(f_scalar2d);
  scalar2d_free(g_scalar2d);
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);
  scalar3d_free(field_scalar3d);

  scalar3d_free(j1_scalar3d);
  scalar3d_free(j2_scalar3d);
  scalar3d_free(j3_scalar3d);
  scalar3d_free(rby_xiplusbbya_scalar3d);
  scalar3d_free(df_dxi_scalar3d);
  scalar3d_free(onebyr_df_dxi_scalar3d);
  scalar3d_free(d2r_dxi2_scalar3d);
  scalar3d_free(onebyr_d2f_dthetadxi_scalar3d);
  scalar3d_free(onebyr_d2r_dthetadxi_scalar3d);
  scalar3d_free(onebyrsintheta_d2f_dphidxi_scalar3d);
  scalar3d_free(onebyrsintheta_d2r_dphidxi_scalar3d);
  scalar3d_free(onebyrsq_anglaplacef_scalar3d);
  scalar3d_free(onebyrsq_anglaplacer_scalar3d);

  scalar3d_free(residual_scalar3d);

  return 0;
}


double ffunc(int z, double theta, double phi)
{
  /* return -3.0 + 1.6*cos(theta); */
  if (z==0)
    return -1.0 - 0.7725484040463793*cos(theta)*sin(theta)*sin(phi);
  else if (z==1)
    return -2.0 - 0.7725484040463793*cos(theta)*sin(theta)*(cos(phi) + sin(phi));
  else
    return 2.0 - 50.0 / (10.0 - 1.5450968080927585*cos(theta)*sin(theta)*(cos(phi) + sin(phi)));
}

double gfunc(int z, double theta, double phi)
{
  if (z==0)
    return 4.0 - 0.7725484040463793*cos(theta)*sin(theta)*cos(phi);
  else if (z==1)
    return 1.0;
  else
    return 0.0;
}

double field(int z, double xi, double theta, double phi)
{
  /*return xi*xi*(3.0*cos(theta)*cos(theta) - 1.0);*/
  return xi*xi*(xi*xi-1)*sin(theta)*cos(theta)*(cos(phi) + sin(phi));
}
