/*********************************************************************
 * remapgridtest.c:                                                  *
 * Tests the functions that remap scalars defined at gridpoints from *
 * one surface matched coordinate system to another:                 *
 * mainly the function "remapgrid" and the functions it calls.       *
 *                                                                   *
 * Author: Benjamin D. Lackey                                        *
 *********************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/poisson.h remapgridtest.c */

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

double field(int z, double xi, double theta, double phi);
void functiontogrid(scalar3d *func_scalar3d, double (*func)(int z, double xi, double theta, double phi));
double fofrtp(int z, double r, double theta, double phi);
void fofrtptogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi));
double T_n(int n, double x);
double boundaryold(int z, double theta, double phi);
double boundarynew(int z, double theta, double phi);

int main (void)
{
  int z, i, j, k;
  int nz = 3;
  int nr = 41; /* must be odd? */
  int nt;
  int np = 40; /* must be even */
  double xi, theta, phi;
  double r_i, theta_j, phi_k;
  double fclenshaw;
  double fdirectsum;
  double ffunc;
  double finterp;
  double error;
  gsl_vector *coeff_vector;
  coeff *f_coeff;
  scalar3d *f_scalar3d;

  scalar2d *boundaryold_scalar2d;
  gsl_vector *alphaold_vector;
  gsl_vector *betaold_vector;
  scalar2d *fold_scalar2d;
  scalar2d *gold_scalar2d;

  scalar2d *boundarynew_scalar2d;
  gsl_vector *alphanew_vector;
  gsl_vector *betanew_vector;
  scalar2d *fnew_scalar2d;
  scalar2d *gnew_scalar2d;
  
  scalar3d *funcnew_scalar3d;
  scalar3d *funcold_scalar3d;
  scalar3d *rnew_scalar3d;

  nt = np/2 + 1;
  
  coeff_vector = gsl_vector_alloc(nr); /* {a_0, a_1, a_2, a_3, a_4, a_5, a_6} */
  f_coeff = coeff_alloc(nz, nr, nt, np);
  f_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  boundaryold_scalar2d = scalar2d_alloc(nz-1, nt, np);
  alphaold_vector = gsl_vector_alloc(nz);
  betaold_vector = gsl_vector_alloc(nz);
  fold_scalar2d = scalar2d_alloc(nz, nt, np);
  gold_scalar2d = scalar2d_alloc(nz, nt, np);

  boundarynew_scalar2d = scalar2d_alloc(nz-1, nt, np);
  alphanew_vector = gsl_vector_alloc(nz);
  betanew_vector = gsl_vector_alloc(nz);
  fnew_scalar2d = scalar2d_alloc(nz, nt, np);
  gnew_scalar2d = scalar2d_alloc(nz, nt, np);

  funcnew_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  funcold_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  rnew_scalar3d = scalar3d_alloc(nz, nr, nt, np);

/*   /\*>>>>>>>>>>>>>>>> test Clenshaw recursion function <<<<<<<<<<<<<<<<*\/ */
  
/*   gsl_vector_set_zero(coeff_vector); */
/*   for ( i = 1; i < nr; i += 2 ) { */
/*     gsl_vector_set(coeff_vector, i, 1.0); */
/*   } */
  
/* /\*   for ( i = 0; i < nr; i++ ) { *\/ */
/* /\*     gsl_vector_set(coeff_vector, i, neg1toi(i)*1.0); *\/ */
/* /\*     /\\*gsl_vector_set(coeff_vector, i, pow(2.0, -1.0*i));*\\/ *\/ */
/* /\*   } *\/ */
/*   print_vector(coeff_vector); */
  
/*   /\* compare Clenshaw recursion relation to a direct summation of Chebyshew polynomials *\/ */
/*   for ( xi = -1.0; xi <= 1.0; xi += 0.1 ) { */
/*     fclenshaw = chebyshevinterpolation(coeff_vector, xi); */
    
/*     fdirectsum = 0.0; */
/*     for ( i = nr-1; i >= 0; i-- ) { */
/*       fdirectsum += gsl_vector_get(coeff_vector, i)*T_n(i, xi); */
/*     } */
/*     /\* relative error *\/ */
/*     printf("(%.18e, %.18e, %.18e, %.18e)\n", xi, fclenshaw, fdirectsum, (fclenshaw-fdirectsum)/fclenshaw); */
/*   } */
  
/*   gsl_vector_free(coeff_vector); */
  
/*   /\*>>>>>>>>>>>>> test spectral interpolation function <<<<<<<<<<<<<<<*\/ */
  
/*   functiontogrid(f_scalar3d, field); */
/*   gridtofourier(f_coeff, f_scalar3d, 0, 0); */
  
/* /\*   coeff_set(f_coeff, 0, 0, 1, 0, 0, 1.0); *\/ */
/* /\*   coeff_set(f_coeff, 1, 1, 1, 0, 0, 1.0); *\/ */
  
/*   print_coeff(f_coeff); */
  
/*   for ( z = 0; z < nz; z++ ){ */
/*     for ( xi = -1.0; xi <= 1.0; xi += 0.2 ) { */
/*       for ( theta = 0; theta <= PI; theta += PI/5.0 ) { */
/* 	for ( phi = 0; phi <= 2.0*PI; phi += 2*PI/5.0 ) { */
/* 	  ffunc = field(z, xi, theta, phi); */
/* 	  finterp = spectralinterpolation(f_coeff, z, xi, theta, phi); */
/* 	  error = finterp - ffunc; */
/* 	  /\*error = (finterp - ffunc) / ffunc;*\/ */
/* 	  /\*error = 2.0*(finterp - ffunc) / (finterp + ffunc);*\/ */
/* 	  printf("z=%d, xi=%.18e, theta=%.18e, phi=%.18e, %.18e, %.18e, %.18e\n", z, xi, theta, phi, finterp, ffunc, error); */
/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   scalar3d_free(f_scalar3d); */
/*   coeff_free(f_coeff); */
  
/*   /\*>>>>>>>>>>>>>>>> test root finder to get xi from r=r(xi) <<<<<<<<<<<<<*\/ */

/*   printf("%.18e\n", xiofroru(3, 1, 1.0, 5.0, 3.0, 2.0, 7.2)); */
/*   for(i = 0; i <= 10; i++ ) { */
/*     xi = i/10.0; */
/*     printf("z=0, %.18e, %.18e\n", xi, xiofroru(3, 0, 1.0, 5.0, 3.0, 2.0, roruofxi(3, 0, 1.0, 5.0, 3.0, 2.0, xi))); */
/*   } */
/*   for(i = 0; i <= 10; i++ ) { */
/*     xi = 2.0*i/10.0 - 1.0; */
/*     printf("z=1, %.18e, %.18e\n", xi, xiofroru(3, 1, 1.0, 5.0, 3.0, 2.0, roruofxi(3, 1, 1.0, 5.0, 3.0, 2.0, xi))); */
/*   } */
/*   for(i = 0; i <= 10; i++ ) { */
/*     xi = 2.0*i/10.0 - 1.0; */
/*     printf("z=2, %.18e, %.18e\n", xi, xiofroru(3, 2, -1.0, 5.0, 1.0, 2.0, roruofxi(3, 2, -1.0, 5.0, 1.0, 2.0, xi))); */
/*   } */
  
  /*>>>>>>>>>>>>>>>>> test remap function <<<<<<<<<<<<<<<<*/
  
  /* Evaluate analytic boundary function on the boundary grid points. */
  for ( z = 0; z < nz-1; z++ ) {
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	scalar2d_set(boundaryold_scalar2d, z, j, k, boundaryold(z, theta_j, phi_k));
	scalar2d_set(boundarynew_scalar2d, z, j, k, boundarynew(z, theta_j, phi_k));
      }
    }
  }
   
  /* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d */
  map_physicaltogrid(boundaryold_scalar2d, alphaold_vector, betaold_vector, fold_scalar2d, gold_scalar2d);
  map_physicaltogrid(boundarynew_scalar2d, alphanew_vector, betanew_vector, fnew_scalar2d, gnew_scalar2d);

  printf("fold_scalar2d\n");
  print_scalar2d(fold_scalar2d);
  printf("gold_scalar2d\n");  
  print_scalar2d(gold_scalar2d);
  printf("alphaold_scalar2d\n");  
  print_vector(alphaold_vector);
  printf("betaold_scalar2d\n");
  print_vector(betaold_vector);
  
  printf("fnew_scalar2d\n");
  print_scalar2d(fnew_scalar2d);
  printf("gnew_scalar2d\n");  
  print_scalar2d(gnew_scalar2d);
  printf("alphanew_scalar2d\n");  
  print_vector(alphanew_vector);
  printf("betanew_scalar2d\n");
  print_vector(betanew_vector);
  
  fofrtptogrid(funcold_scalar3d, alphaold_vector, betaold_vector, fold_scalar2d, gold_scalar2d, fofrtp);

  printf("hi there\n");
  remapgrid(funcnew_scalar3d, alphanew_vector, betanew_vector, fnew_scalar2d, gnew_scalar2d, funcold_scalar3d, alphaold_vector, betaold_vector, fold_scalar2d, gold_scalar2d);
  printf("hi there 2\n");

  rofxtp(rnew_scalar3d, alphanew_vector, betanew_vector, fnew_scalar2d, gnew_scalar2d);
  
   /* print to file the position and function value at all gridpoints */
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  r_i = scalar3d_get(rnew_scalar3d, z, i, j, k);
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  ffunc = fofrtp(z, r_i, theta_j, phi_k);
	  finterp = scalar3d_get(funcnew_scalar3d, z, i, j, k);
	  /*error = finterp - ffunc;*/
	  error = (finterp - ffunc) / ffunc;
	  /*error = 2.0*(finterp - ffunc) / (finterp + ffunc);*/
	  printf("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", r_i, theta_j, phi_k, ffunc, finterp, error);
	}
      }
    }
  }
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	r_i = 1.0 / scalar3d_get(rnew_scalar3d, z, i, j, k);
	theta_j = PI*j/(nt-1);
	phi_k = 2*PI*k/np;
	ffunc = fofrtp(z, r_i, theta_j, phi_k);
	finterp = scalar3d_get(funcnew_scalar3d, z, i, j, k);
	/*error = finterp - ffunc;*/
	error = (finterp - ffunc) / ffunc;
	/*error = 2.0*(finterp - ffunc) / (finterp + ffunc);*/
	printf("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", r_i, theta_j, phi_k, ffunc, finterp, error);
      }
    }
  }
    

  /*scalar2d_free(boundaryold_scalar2d);*/
  gsl_vector_free(alphaold_vector);
  gsl_vector_free(betaold_vector);
  /*scalar2d_free(fold_scalar2d);*/
  /*scalar2d_free(gold_scalar2d);*/
  
  /*scalar2d_free(boundarynew_scalar2d);*/
  gsl_vector_free(alphanew_vector);
  gsl_vector_free(betanew_vector);
  /*scalar2d_free(fnew_scalar2d);*/
  /*scalar2d_free(gnew_scalar2d);*/
  
  scalar3d_free(funcnew_scalar3d);
  scalar3d_free(funcold_scalar3d);
  scalar3d_free(rnew_scalar3d);
  
  return 0;
}


/* double boundaryold(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.5; */
/*   else */
/*     return 3.0; */
/* } */

double boundaryold(int z, double theta, double phi)
{
  if(z==0)
    return 2.0*(1.0 + /*0.3*sin(theta)*(cos(phi)+sin(phi)) +*/ 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi)));
  else
    return 5.0;
}

/* double boundaryold(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else */
/*     return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/* } */

double boundarynew(int z, double theta, double phi)
{
  if(z==0)
    return 1.0;
  else
    return 5.0;
}

/* double field(int z, double xi, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return sin(xi)*sin(sin(theta))*(exp(cos(phi))+exp(sin(phi)))  */
/*       + cos(xi)*cos(cos(theta))*(exp(cos(phi))+exp(sin(phi))); */
/*   else */
/*     return exp(xi)*(exp(sin(theta))*(exp(cos(phi))+exp(sin(phi))) + cos(cos(theta))*(exp(cos(phi))+exp(sin(phi)))); */
/* } */


double field(int z, double xi, double theta, double phi)
{
  if(z==0)
    return 2.222*T_n(2, xi)*cos(6*theta)*cos(6*phi);
  else
    return 2.222*T_n(1, xi)*cos(6*theta)*cos(6*phi);
}

/* double fofrtp(int z, double r, double theta, double phi) */
/* { */
/*   if(z<2) */
/*     return r*r*exp(-r*r/5.0); */
/*   else */
/*     return 1/(r*r); */
/* } */

double fofrtp(int z, double r, double theta, double phi)
{
  return r*r*exp(-r*r/5.0);
}

/***************************************************************************************/
/* Take a function f = f(z, xi, theta, phi) and evaluate it on the grid func_scalar3d. */
/***************************************************************************************/ 
void functiontogrid(scalar3d *func_scalar3d, double (*func)(int z, double xi, double theta, double phi))
{
  int z, i, j, k;
  int nz, nr, nt, np;
  double xi_i, theta_j, phi_k;
  
  nz = func_scalar3d->nz;
  nr = func_scalar3d->nr;
  nt = func_scalar3d->nt;
  np = func_scalar3d->np;
  
  /* kernel and shells */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if(z==0) /* kernel */
	    xi_i = sin(PI*i/(2*(nr-1)));
	  else /* shells or external domain */
	    xi_i = -cos(PI*i/(nr-1));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  scalar3d_set(func_scalar3d, z, i, j, k, func(z, xi_i, theta_j, phi_k));
	}
      }
    }
  }
}


/******************************************************************************************************/
/* Take a function f = f(z, r, theta, phi) and evaluate it on the surface matched grid func_scalar3d. */
/******************************************************************************************************/ 
void fofrtptogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi))
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


/************************************************/
/* Evaluate nth Chebyshev polynomial at point x */
/************************************************/
double T_n(int n, double x)
{
  double Tn;          /* nth Chebyshev polynomial T_n(x) */
  double Tnminus1;    /* (n-1)th Chebyshev polynomial T_{n-1}(x) */ 
  double Tnminus2;    /* (n-2)th Chebyshev polynomial T_{n-2}(x) */ 
  int j; 
    
  /* Evaluate nth Chebyshev polynomial at point x.     */
  /* Use the recursion relation:                       */
  /* T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x)              */
  if(n<0) {
    printf("can't have n<0 in T_n\n");
    return 0.0;
  } else if(n==0) {
    return 1;
  } else if(n==1) {
    return x;
  } else {
    Tnminus2 = 1;
    Tnminus1 = x;
    Tn = 0;
    for(j=2; j<=n; j++) {
      Tn = 2.0*x*Tnminus1 - Tnminus2;
      Tnminus2 = Tnminus1;
      Tnminus1 = Tn;
    }
    return Tn;
  } 
}
