/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/poisson.h j2test.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

double boundary(int z, int nt, int np, double theta, double phi);
double dbounddt(int z, int nt, int np, double theta, double phi);
double field(int z, int nr, int nt, int np, double xi, double theta, double phi);
double fieldphysical(int z, double r, double theta, double phi);
double T_n(int n, double x);

int main (void)
{

  int z, j, k, imag;
  int nz = 3;
  int nt = 15;
  int np = 28; /* must be even */
  int npc;
  double xi_i, theta_j, phi_k;
  scalar2d *f_grid;
  bound_coeff *f_fourier;
  scalar2d *dfdt_grid;
  scalar2d *dfdtanalytic_grid;
  
  npc = ( np / 2 ) + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  f_grid = scalar2d_alloc(nz, nt, np);
  f_fourier = bound_coeff_alloc(nz, nt, np);
  dfdt_grid = scalar2d_alloc(nz, nt, np);
  dfdtanalytic_grid = scalar2d_alloc(nz, nt, np);
  
  /* Evaluate analytic 2d function on the grid points. */
  for ( z = 0; z < nz; z++ ) {
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	scalar2d_set(f_grid, z, j, k, boundary(z, nt, np, theta_j, phi_k));
      }
    }
  }
  print_scalar2d(f_grid);
  
  gridtofourier_bound(f_fourier, f_grid);
  print_bound_coeff(f_fourier);
  
  /* d/dt[a_j cos(j theta)] = -j a_j sin(j theta) */  
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) {   
      for(k=2*imag; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j=0; j<nt; j++) {
	  /* f: */
	  bound_coeff_set(f_fourier, z, j, k, imag, 
			  -j * bound_coeff_get(f_fourier, z, j, k, imag));
	}
      }
    }
  }
  
  /* d/dt[a_j sin(j theta)] = +j a_j cos(j theta) */
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) {   
      for(k=1; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j=1; j<nt-1; j++) {  
	  /* f: */
	  bound_coeff_set(f_fourier, z, j, k, imag, 
			  j * bound_coeff_get(f_fourier, z, j, k, imag));
	}
      }
    }
  }
  print_bound_coeff(f_fourier);
  
  fouriertogrid_bound(dfdt_grid, f_fourier, 1, 0);
  print_scalar2d(dfdt_grid);
  
  /* Evaluate analytic 2d function on the grid points. */
  for ( z = 0; z < nz; z++ ) {
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	scalar2d_set(dfdtanalytic_grid, z, j, k, dbounddt(z, nt, np, theta_j, phi_k));
      }
    }
  }
  print_scalar2d(dfdtanalytic_grid);
  
  /* subtract analytical from numerical derivative. */
  for ( z = 0; z < nz; z++ ) {
    for ( k = 0; k < np; k++ ) {
      for ( j = 0; j < nt; j++ ) {
	scalar2d_set(dfdtanalytic_grid, z, j, k, scalar2d_get(dfdt_grid, z, j, k) - scalar2d_get(dfdtanalytic_grid, z, j, k));
      }
    }
  }
  print_scalar2d(dfdtanalytic_grid);

  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0; */
/*   else */
/*     return 5.0; */
/* } */

/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else if(z==1) */
/*     return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/* } */

double boundary(int z, int nt, int np, double theta, double phi)
{
  int j, k;
  double sum;
  int npc = np/2 + 1;

/*   sum = 0; */
/*   for ( j = 0; j < nt-1; j++ ) { */
/*     for ( k = 0; k < npc; k+=2 ) { */
/*       sum += 1.1111111*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/*     } */
/*     for ( k = 1; k < npc; k+=2 ) { */
/*       sum += 1.1111111*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/*     } */
/*   } */
  sum = 0;
  for ( j = 0; j < 5; j++ ) {
    for ( k = 0; k < 5; k+=2 ) {
      sum += 1.1111111*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
    }
    for ( k = 1; k < 5; k+=2 ) {
      sum += 1.1111111*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
    }
  }
  return sum;
}

double dbounddt(int z, int nt, int np, double theta, double phi)
{
  int j, k;
  double sum;
  int npc = np/2 + 1;

/*   sum = 0; */
/*   for ( j = 0; j < nt-1; j++ ) { */
/*     for ( k = 0; k < npc; k+=2 ) { */
/*       sum += 1.1111111*(-j)*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/*     } */
/*     for ( k = 1; k < npc; k+=2 ) { */
/*       sum += 1.1111111*(j)*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/*     } */
/*   } */
  sum = 0;
  for ( j = 0; j < 5; j++ ) {
    for ( k = 0; k < 5; k+=2 ) {
      sum += 1.1111111*(-j)*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
    }
    for ( k = 1; k < 5; k+=2 ) {
      sum += 1.1111111*(j)*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
    }
  }
  return sum;
}


/*********************************************************/
/* Some function of position in spherical coordinates.   */
/*********************************************************/
double fieldphysical(int z, double r, double theta, double phi)
{
  return cos(r);
}


/******************************************************/
/* Evaluate a function at the point (xi, theta, phi). */
/******************************************************/
double field(int z, int nr, int nt, int np, double xi, double theta, double phi)
{
  return cos(xi);
}

/* double field(int z, int nr, int nt, int np, double xi, double theta, double phi) */
/* { */
/*   int i; */
/*   double sum; */
  
/*   if(z==0) { */
/*     sum = 0.0; */
/*     for(i=0; i<nr; i++){ */
/*       sum += T_n(2*i, xi); */
/*     } */
/*     for(i=0; i<nr-1; i++){ */
/*       sum += T_n(2*i+1, xi)*(cos(theta) + sin(theta)*(cos(phi)+sin(phi))); */
/*     } */
/*     return sum; */
/*   } else { */
/*     sum = 0.0; */
/*     for(i=0; i<nr; i++){ */
/*       sum += T_n(i, xi); */
/*     } */
/*     return sum; */
/*   } */
/* } */


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

  if(n==0) {
	return 1;
  } else if(n==1) {
	return x;
  } else {
	Tnminus2 = 1;
	Tnminus1 = x;
	for(j=2; j<=n; j++) {
	  Tn = 2.0*x*Tnminus1 - Tnminus2;
	  Tnminus2 = Tnminus1;
	  Tnminus1 = Tn;
	}
	return Tn;
  } 
}
