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

void randomboundary(scalar2d *boundary_scalar2d);
void errorlist(scalar3d *f_num_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi));
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
  
  double grad, grad_analytic, error, max;

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
  /*randomboundary(boundary_scalar2d);*/
  print_scalar2d(boundary_scalar2d);
  
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
  max = 0.0;
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
	  if((z == nz-1)&&(i<nr-1)&&(ABS(error)<1.0e-3))
	    max = MAX(error, max);
	}
      }
    }
  }
  printf("Max error in external zone = %.18e\n", max);

  /*errorlist(gradf_phi_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, gradfield_phi);*/

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


void randomboundary(scalar2d *boundary_scalar2d)
{
  int z, j, k;
  int nz, nt, np;
  double bound;
  
  nz = boundary_scalar2d->nz; /* Here, nz is number of boundaries not zones. */
  nt = boundary_scalar2d->nt;
  np = boundary_scalar2d->np;

  srand(3); /* seed for random number generator rand() */

  for ( z = 0; z < nz; z++ ) {
    /* all points on N pole must be same */
    j = 0;
    bound = 5.0*(z+1) + 0.3*((double)rand()/(double)RAND_MAX - 0.5);
    for ( k = 0; k < np; k++ )
      scalar2d_set(boundary_scalar2d, z, j, k, bound);
    
    for ( j = 1; j < nt-1; j++ ) {
      for ( k = 0; k < np; k++ ) {
	bound = 5.0*(z+1) + 0.3*((double)rand()/(double)RAND_MAX - 0.5);
	scalar2d_set(boundary_scalar2d, z, j, k, bound);
      }
    }
    
    j = nt-1;
    bound = 5.0*(z+1) + 0.3*((double)rand()/(double)RAND_MAX - 0.5);
    for ( k = 0; k < np; k++ )
      scalar2d_set(boundary_scalar2d, z, j, k, bound);
  }
}


void errorlist(scalar3d *f_num_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi))
{
  int z, i, j, k;
  int nz, nr, nt, np;
  scalar3d *r_scalar3d;
  double r_i, theta_j, phi_k;
  double f_anal, f_num, f_max;
  int p;
  double error, pnorm;
  
  nz = f_num_scalar3d->nz;
  nr = f_num_scalar3d->nr;
  nt = f_num_scalar3d->nt;
  np = f_num_scalar3d->np;
  
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* Find max value of f */
  f_max = 0.0;
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  f_anal = func(z, r_i, theta_j, phi_k);
	  f_max = MAX(f_anal, f_max);
	}
      }
    }
  }

  /* List the error and find the p-norm */
  p = 1;
  pnorm = 0.0;
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	printf("-----------------------------------------\n");
	for ( k = 0; k < np; k++ ) {
	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  f_num = scalar3d_get(f_num_scalar3d, z, i, j, k);
	  f_anal = func(z, r_i, theta_j, phi_k);
	  error = MIN( ABS((f_num - f_anal)/f_anal), ABS(f_num/f_max) );
	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, r_i, theta_j, phi_k, f_num, f_anal, error);
	  pnorm += pow(ABS(error), p); 
	}
      }
    }
  }
  pnorm = pow(pnorm, 1.0/p) / (nz*nr*nt*np);
  printf("%d-norm = %.18e\n", p, pnorm);

  scalar3d_free(r_scalar3d);
}

/******************************/
/* Function for the boundary. */
/******************************/
/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0; */
/*   else */
/*     return 5.0; */
/* } */
double boundary(int z, double theta, double phi)
{
  if(z==0)
    return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/* + 0.2*cos(5*theta)*/);
  else
    return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))/* + 0.2*cos(5*theta)*/);
}


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
    return r*r*sin(theta)*sin(theta)*(cos(2.0*phi) + sin(2.0*phi));
  else
    return sin(theta)*sin(theta)*(cos(2.0*phi) + sin(2.0*phi))/(r*r);
}

/**************************************/
/* The corresponding phi gradient.    */
/**************************************/
double gradfield_phi(int z, double r, double theta, double phi)
{
  if(z<2)
    return r*sin(theta)*(-2.0*sin(2.0*phi) + 2.0*cos(2.0*phi));
  else
    return sin(theta)*(-2.0*sin(2.0*phi) + 2.0*cos(2.0*phi))/(r*r*r);
}
