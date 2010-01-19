/**************************************************************
 * gradienttest.c:                                            *
 * Tests functions that calculate the gradient in             *
 * spherical coordinates at surface-matched gridpoints.       *
 *                                                            *
 * gradient_r tests the functions:                            *
 *   gridtofourier                                            *
 *   dfdxi                                                    *
 *   fouriertogrid                                            *
 *   jacobian1                                                *
 *   rofxtp                                                   *
 *                                                            *
 * gradient_theta tests the additional functions:             *
 *   dfdthetaprime                                            *
 *   dividebyr                                                *
 *   jacobian2                                                *
 *                                                            *
 * gradient_phi tests the additional functions:               *
 *   dfdphiprime                                              *
 *   dividebysin                                              *
 *   dividebyr                                                *
 *   jacobian3                                                *
 *                                                            *
 * Author: Benjamin D. Lackey                                 *
 **************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/poisson.h gradienttest.c */

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
double fieldforr(int z, double r, double theta, double phi);
double gradfield_r(int z, double r, double theta, double phi);
double fieldfortheta(int z, double r, double theta, double phi);
double gradfield_theta(int z, double r, double theta, double phi);
double fieldforphi(int z, double r, double theta, double phi);
double gradfield_phi(int z, double r, double theta, double phi);


int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 35; /* must be odd? */
  int nt;
  int np = 30; /* must be even */
  double r_i, theta_j, phi_k;
  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *field_scalar3d;
  scalar3d *gradf_r_scalar3d; 
  scalar3d *gradf_theta_scalar3d;
  scalar3d *gradf_phi_scalar3d;
  scalar3d *r_scalar3d;
  
  double grad, grad_analytic, error;

  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  gradf_r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  gradf_theta_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  gradf_phi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* evaluate boundary function on gridpoints */
  boundarytogrid(boundary_scalar2d, boundary);
  /*print_scalar2d(boundary_scalar2d);*/

  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
/*   print_vector(alpha_vector); */
/*   print_vector(beta_vector); */
/*   print_scalar2d(f_scalar2d); */
/*   print_scalar2d(g_scalar2d); */

  /*>>>>>>>>>>>>>>>>>> TESTING GRAD_R <<<<<<<<<<<<<<<<<<*/

/*   /\* evaluate field at surface matched gridpoints *\/ */
/*   functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, fieldforr); */

/*   /\* calculate radial gradient *\/ */
/*   gradient_r(gradf_r_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */

/*   /\* Find the radial position of each point. *\/ */
/*   rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */
/*   print_scalar3d(r_scalar3d); */
/*   /\* compare numerical to analytic values of the gradient at the surface matched gridpoints *\/ */
/*   printf("GRAD_R AT SURFACE MATCHED GRIDPOINTS:\n"); */
/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k)); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/* 	  grad = scalar3d_get(gradf_r_scalar3d, z, i, j, k); */
/* 	  grad_analytic = gradfield_r(z, r_i, theta_j, phi_k); */
/* 	  error = (grad - grad_analytic)/grad_analytic; */
/* 	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, r_i, theta_j, phi_k, grad, grad_analytic, error); */
/* 	} */
/*       } */
/*     } */
/*   } */
  
  
  /*>>>>>>>>>>>>>>>>>> TESTING GRAD_THETA <<<<<<<<<<<<<<<<<<*/
  
/*   /\* evaluate field at surface matched gridpoints *\/ */
/*   functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, fieldfortheta); */
/* /\*   print_scalar3d(field_scalar3d); *\/ */

/*   /\* calculate radial gradient *\/ */
/*   gradient_theta(gradf_theta_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */
  
/*   /\* Find the radial position of each point. *\/ */
/*   rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */
/*   /\*print_scalar3d(r_scalar3d);*\/ */
/*   /\* compare numerical to analytic values of the gradient at the surface matched gridpoints *\/ */
/*   printf("GRAD_THETA AT SURFACE MATCHED GRIDPOINTS:\n"); */
/*   for ( z = 0; z < nz; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k)); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/* 	  grad = scalar3d_get(gradf_theta_scalar3d, z, i, j, k); */
/* 	  grad_analytic = gradfield_theta(z, r_i, theta_j, phi_k); */
/* 	  error = (grad - grad_analytic)/(grad); */
/* 	  /\*error = MIN(ABS((grad - grad_analytic)/grad_analytic), ABS((grad - grad_analytic)/(grad + grad_analytic+1.0e-15)));*\/ */
/* 	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, r_i, theta_j, phi_k, grad, grad_analytic, error); */
/* 	} */
/*       } */
/*     } */
/*   } */
 

  /*>>>>>>>>>>>>>>>>>> TESTING GRAD_PHI <<<<<<<<<<<<<<<<<<*/
  
  /* evaluate field at surface matched gridpoints */
  functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, fieldforphi);
  
  /* calculate radial gradient */
  gradient_phi(gradf_phi_scalar3d, field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* compare numerical to analytic values of the gradient at the surface matched gridpoints */
  printf("GRAD_PHI AT SURFACE MATCHED GRIDPOINTS:\n");
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	printf("-----------------------------------------\n");
	for ( k = 0; k < np; k++ ) {
	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  grad = scalar3d_get(gradf_phi_scalar3d, z, i, j, k);
	  grad_analytic = gradfield_phi(z, r_i, theta_j, phi_k);
	  error = (grad - grad_analytic)/(grad);
	  /*error = MIN(ABS((grad - grad_analytic)/grad_analytic), ABS((grad - grad_analytic)/(grad + grad_analytic+1.0e-15)));*/
	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, r_i, theta_j, phi_k, grad, grad_analytic, error);
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
  scalar3d_free(gradf_r_scalar3d);
  scalar3d_free(gradf_theta_scalar3d);
  scalar3d_free(gradf_phi_scalar3d);
  scalar3d_free(r_scalar3d);

  return 0;
}

/******************************/
/* Function for the boundary. */
/******************************/
double boundary(int z, double theta, double phi)
{
  if(z==0)
    return 1.0;
  else
    return 5.0;
}
/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else */
/*     return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/* } */


/*>>>>>>>>>>>>> RADIAL <<<<<<<<<<<<<<*/

/*********************************************/
/* Some function in spherical coordinates.   */
/*********************************************/
double fieldforr(int z, double r, double theta, double phi)
{
  if(z<2)
    return cos(0.2*r);
  else
    return 1.0/(r*r);
}

/**************************************/
/* The corresponding radial gradient. */
/**************************************/
double gradfield_r(int z, double r, double theta, double phi)
{
  if(z<2)
    return -0.2*sin(0.2*r);
  else
    return -2.0/(r*r*r);
}

/*>>>>>>>>>>>>> THETA <<<<<<<<<<<<<<*/

/*********************************************/
/* Some function in spherical coordinates.   */
/*********************************************/
double fieldfortheta(int z, double r, double theta, double phi)
{
  if(z<2)
    return r*cos(theta);
  else
    return cos(theta)/r;
}

/**************************************/
/* The corresponding theta gradient.  */
/**************************************/
double gradfield_theta(int z, double r, double theta, double phi)
{
  if(z<2)
    return -sin(theta);
  else
    return -sin(theta)/(r*r);
}

/*>>>>>>>>>>>>> PHI <<<<<<<<<<<<<<*/

/*********************************************/
/* Some function in spherical coordinates.   */
/***** ***************************************/
double fieldforphi(int z, double r, double theta, double phi)
{
  if(z<2)
    return r*sin(theta)*(cos(phi) + sin(phi));
  else
    return sin(theta)*(cos(phi) + sin(phi))/r;
}

/**************************************/
/* The corresponding phi gradient.    */
/**************************************/
double gradfield_phi(int z, double r, double theta, double phi)
{
  if(z<2)
    return -sin(phi) + cos(phi);
  else
    return (-sin(phi) + cos(phi))/(r*r);
}
