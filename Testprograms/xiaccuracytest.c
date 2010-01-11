/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/poisson.h xiaccuracytest.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

/*double boundary(int z, double theta, double phi); */
/* double boundary(int z, int nt, int np, double theta, double phi); */
/* double field(int z, int nr, int nt, int np, double xi, double theta, double phi); */
double field(int z, int nr, double xi);
/*double fieldphysical(double r, double theta, double phi);*/
double T_n(int n, double x);

int main (void)
{

 
  /*clock_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12; 
    float ratio;*/ 

  FILE *fpgrid;
  /*FILE *fpfunction;*/

  int z, i, j, k;
  int nz = 5;
  int nr = 5; /* must be odd? */
  int nt = 7;
  int np = 12; /* must be even */
  int npc;
  double xi_i;
  scalar3d *field_grid;
  coeff *field_coeff;
  coeff *dfdxi_coeff;
  scalar3d *dfdxi_grid;
  scalar3d *temp_grid;
  double f_zijk;
  double diff;

  npc = ( np / 2 ) + 1;
  
  /* the field */
  field_grid = scalar3d_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);
  
  /* derivatives of field */
  dfdxi_grid = scalar3d_alloc(nz, nr, nt, np);
  dfdxi_coeff = coeff_alloc(nz, nr, nt, np); 
  
  temp_grid = scalar3d_alloc(nz, nr, nt, np); 

  /* Assign the field data to field_grid. */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if(z==0)
	    xi_i = sin(PI*i/(2*(nr-1)));
	  else
	    xi_i = -cos(PI*i/(nr-1));
	  scalar3d_set(field_grid, z, i, j, k, field(z, nr, xi_i));
	}
      }
    }
  }
    
  gridtofourier(field_coeff, field_grid, 0, 0); 
  dfdxi(dfdxi_coeff, field_coeff);
  fouriertogrid(dfdxi_grid, dfdxi_coeff, 1, 0);
 
  /* print to file the function and radial position for all gridpoints */
  fpgrid=fopen("fofr.txt", "w");
  for ( z = 1; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	   if(z==0)
	    xi_i = sin(PI*i/(2*(nr-1)));
	  else
	    xi_i = -cos(PI*i/(nr-1));
	  f_zijk = scalar3d_get(dfdxi_grid, z, i, j, k);
	  diff = (scalar3d_get(dfdxi_grid, z, i, j, k) - (1.0))/(1.0);
	  fprintf(fpgrid, "%.18e\t%.18e\t%.18e\n", xi_i, f_zijk , diff);
	}
      }
    }
  }
  
  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
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
/*     return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else */
/*     return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/* } */

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
/* double fieldphysical(double r, double theta, double phi) */
/* { */
/*   return cos(0.2*r); */
/* } */


/* sum = 1.0; */
/* for(i=0; i<nr; i++){ */
/*   newsum = xi*(1.0+sum); */
/*   sum = newsum; */
/*  } */
/******************************************************/
/* Evaluate a function at the point (xi, theta, phi). */
/******************************************************/
double field(int z, int nr, double xi)
{
  return xi;
  /*  return pow(xi, 6);*/
/*   int i; */
/*   int sum; */
  
/*   if(z==0) { */
/*     sum = 0.0; */
/*     for(i=0; i<nr; i+=2) { */
/*       sum += pow(xi, i); */
/*     } */
/*     return sum; */
/*   } else { */
/*     sum = 0.0; */
/*     for(i=0; i<nr; i++) { */
/*       sum += pow(xi, i); */
/*     } */
/*     return sum; */
/*   } */
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
