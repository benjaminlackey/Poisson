/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/poisson.h dividebysintest.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

double boundary(int z, int nt, int np, double theta, double phi);
double field(int z, int nr, int nt, int np, double xi, double theta, double phi);
double fieldphysical(int z, double r, double theta, double phi);
double T_n(int n, double x);

int main (void)
{
  int z, i, j, k;
  int imag;
  int nz = 3;
  int nr = 5; /* must be odd? */
  int nt = 7;
  int np = 12; /* must be even */
  int npc;
  coeff *f_coeff;
  coeff *fbysin_coeff;
  bound_coeff *f_bound_coeff;
  bound_coeff *fbysin_bound_coeff;

  npc = ( np / 2 ) + 1;

  f_coeff = coeff_alloc(nz, nr, nt, np);
  fbysin_coeff = coeff_alloc(nz, nr, nt, np);

  f_bound_coeff = bound_coeff_alloc(nz, nt, np);
  fbysin_bound_coeff = bound_coeff_alloc(nz, nt, np);
 
  /* Fill up coefficient structure. */
  for(imag=0; imag<=1; imag++) {
    for ( z = 0; z < nz; z++ ) {
      for ( i = 0; i < nr; i++ ) {
	/* fill up sin(theta) coefficients */
	for ( k = 1; k < npc-1; k+=2 ) {
	  for ( j = 1; j < nt-1; j++ ) {
	    coeff_set(f_coeff, z, i, j, k, imag, 1.0);
	  }
	}
	/* fill up cos(theta) coefficients*/ 
	/* (currently using coefficients of P^6_2(cos(theta)) + P^5_4(cos(theta))) */
	for ( k = 0; k < npc; k+=2 ) {
	  coeff_set(f_coeff, z, i, 0, k, imag, 10.0*105.0/256.0);
	  coeff_set(f_coeff, z, i, 1, k, imag, 288.0*105.0/256.0);
	  coeff_set(f_coeff, z, i, 2, k, imag, 17.0*105.0/256.0);
	  coeff_set(f_coeff, z, i, 3, k, imag, -432.0*105.0/256.0);
	  coeff_set(f_coeff, z, i, 4, k, imag, 6.0*105.0/256.0);
	  coeff_set(f_coeff, z, i, 5, k, imag, 144.0*105.0/256.0);
	  coeff_set(f_coeff, z, i, 6, k, imag, -33.0*105.0/256.0);
	  
	}
      }
    }
  }
  print_coeff(f_coeff);
  
  dividebysin(fbysin_coeff, f_coeff);
  /*print_coeff(fbysin_coeff);*/


  /*>>>>>>>>>>> now do the same thing for the boundary coefficients <<<<<<<<<<<<<*/
  /* Fill up coefficient structure. */
  for(imag=0; imag<=1; imag++) {
    for ( z = 0; z < nz; z++ ) {
      /* fill up sin(theta) coefficients */
      for ( k = 1; k < npc-1; k+=2 ) {
	for ( j = 1; j < nt-1; j++ ) {
	  bound_coeff_set(f_bound_coeff, z, j, k, imag, 1.0);
	}
      }
      /* fill up sin(theta) coefficients (currently using */
      for ( k = 0; k < npc; k+=2 ) {
	bound_coeff_set(f_bound_coeff, z, 0, k, imag, 10.0*105.0/256.0);
	bound_coeff_set(f_bound_coeff, z, 1, k, imag, 288.0*105.0/256.0);
	bound_coeff_set(f_bound_coeff, z, 2, k, imag, 17.0*105.0/256.0);
	bound_coeff_set(f_bound_coeff, z, 3, k, imag, -432.0*105.0/256.0);
	bound_coeff_set(f_bound_coeff, z, 4, k, imag, 6.0*105.0/256.0);
	bound_coeff_set(f_bound_coeff, z, 5, k, imag, 144.0*105.0/256.0);
	bound_coeff_set(f_bound_coeff, z, 6, k, imag, -33.0*105.0/256.0);
  
      }
    }
  }
  print_bound_coeff(f_bound_coeff);
  
  dividebysin_bound(fbysin_bound_coeff, f_bound_coeff);
  print_bound_coeff(fbysin_bound_coeff);
  
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

double boundary(int nt, int np, int z, double theta, double phi)
{
  if(z==0)
    return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi)));
  else if(z==1)
    return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi)));
}

/* double boundary(int z, int nt, int np, double theta, double phi) */
/* { */
/*   int j, k; */
/*   double sum; */
/*   int npc = np/2 + 1; */

/*   sum = 1; */
/*   for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < npc; k+=2 ) { */
/* 	  sum += /\*((j+1)*10+(k+1))*\/1.1111111*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/* 	} */
/* 	for ( k = 1; k < npc; k+=2 ) { */
/* 	  sum += /\*((j+1)*10+(k+1))*\/1.1111111*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/* 	} */
/*   } */
  
/*   /\*sum = 0; */
/*   for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np/2+1; k+=2 ) { */
/* 	  sum += (j+1)*cos(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	} */
/* 	for ( k = 1; k < np/2+1; k+=2 ) { */
/* 	  sum += (j+1)*sin(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	} */
/* 	}*\/ */
  
/*   /\*sum = 0; */
/* 	for ( j = 0; j < nt; j++ ) { */
/* 	sum += (j+1)*cos(j*theta); */
/* 	}*\/ */
  
/*   return sum; */
/* } */


/*********************************************************/
/* Some function of position in spherical coordinates.   */
/*********************************************************/
double fieldphysical(int z, double r, double theta, double phi)
{
  return r*cos(theta);
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
