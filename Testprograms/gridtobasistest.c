/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/poisson.h gridtobasistest.c */


/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

double boundary(int z, int nt, int np, double theta, double phi);
double field(int z, int nr, int nt, int np, double xi, double theta, double phi);
/*double fieldphysical(int z, double r, double theta, double phi);*/
double T_n(int n, double x);

int main (void)
{
  /*FILE *fpgrid;
    FILE *fpfunction;*/

  int z, i, j, k;
  double xi_i, theta_j, phi_k;
  int nz = 3;
  int nr = 5; /* must be odd? */
  int nt = 7;
  int np = 6; /* must be even */
  int npc;
  scalar2d *b_zjk;
  bound_coeff *bcoeff_zjk;
  scalar3d *c_zijk;
  coeff *coeff_zijk;

  b_zjk = scalar2d_alloc(nz-1, nt, np);
  bcoeff_zjk = bound_coeff_alloc(nz-1, nt, np);

  c_zijk = scalar3d_alloc(nz, nr, nt, np);
  coeff_zijk = coeff_alloc(nz, nr, nt, np);
 

  npc = ( np / 2 ) + 1;
  
  /* Assign the boundary data. */
  for ( z = 0; z < nz-1; z++ ) {
	for ( k = 0; k < np; k++ ) {
	  phi_k = 2*PI*k/np;
	  for ( j = 0; j < nt; j++ ) {
		theta_j = PI*j/(nt-1);
		scalar2d_set(b_zjk, z, j, k, boundary(z, nt, np, theta_j, phi_k));
	  }
	}
  }
  
/*   /\* Assign the boundary data version 2. *\/ */
/*   for ( z = 0; z < nz-1; z++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  scalar2d_set(b_zjk, z, 0, k, (z+1)*100); */
/* 	  for ( j = 1; j < nt-1; j++ ) { */
/* 		scalar2d_set(b_zjk, z, j, k, (z+1)*100 + (j+1)*10 + (k+1)); */
/* 	  } */
/* 	  scalar2d_set(b_zjk, z, nt-1, k, (z+1)*100); */
/* 	} */
/*   } */
  
  /* Assign the field data. */
  for ( z = 0; z < nz; z++ ) {
	for ( i = 0; i < nr; i++ ) {
	  for ( j = 0; j < nt; j++ ) {
		for ( k = 0; k < np; k++ ) {
		  if(z==0) {
			xi_i = sin(PI*i/(2*(nr-1)));
		  } else {
			xi_i = -cos(PI*i/(nr-1));
		  }
		  theta_j = PI*j/(nt-1);
		  phi_k = 2*PI*k/np;
		  scalar3d_set(c_zijk, z, i, j, k, field(z, nr, nt, np, xi_i, theta_j, phi_k));
		}
	  }
	}
  }
  
  /* print start data*/
  /*printf("  Boundary data.\n");*/
  
  /*print_scalar2d(b_zjk);*/
  
  gridtofourier_bound(bcoeff_zjk, b_zjk);
  
  /*print_bound_coeff(bcoeff_zjk);*/
  
  fouriertogrid_bound(b_zjk, bcoeff_zjk);

  /*print_scalar2d(b_zjk);*/

  /*printf("  Field data.\n");*/

  /*print_scalar3d(c_zijk);*/

  gridtofourier(coeff_zijk, c_zijk);

  print_coeff(coeff_zijk);

  fouriertogrid(c_zijk, coeff_zijk);
  
  print_scalar3d(c_zijk);
  
  gridtofourier(coeff_zijk, c_zijk);
  
  print_coeff(coeff_zijk);
  
  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.1*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else */
/*     return 5.0; */
/* } */
double boundary(int z, int nt, int np, double theta, double phi)
{
  int j, k;
  double sum;
  int npc = np/2 + 1;

  sum = 1;
  for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < npc; k+=2 ) {
	  sum += /*((j+1)*10+(k+1))*/1.1111111*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
	}
	for ( k = 1; k < npc; k+=2 ) {
	  sum += /*((j+1)*10+(k+1))*/1.1111111*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
	}
  }
  
  /*sum = 0;
  for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np/2+1; k+=2 ) {
	  sum += (j+1)*cos(j*theta)*(cos(k*phi) + sin(k*phi));
	}
	for ( k = 1; k < np/2+1; k+=2 ) {
	  sum += (j+1)*sin(j*theta)*(cos(k*phi) + sin(k*phi));
	}
	}*/
  
  /*sum = 0;
	for ( j = 0; j < nt; j++ ) {
	sum += (j+1)*cos(j*theta);
	}*/
  
  return sum;
}


/*********************************************************/
/* Some function of position in spherical coordinates.   */
/*********************************************************/
/* double fieldphysical(int z, double r, double theta, double phi) */
/* { */
/*   return cos(r); */
/* } */


/******************************************************/
/* Evaluate a function at the point (xi, theta, phi). */
/******************************************************/
double field(int z, int nr, int nt, int np, double xi, double theta, double phi)
{
  int i, j, k;
  double basis; /* basis function */
  double sum;
  int npc = np/2 + 1;
  
  sum = 0.0;
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < npc; k++ ) {
      
      if( z==0 && j%2 ) { /* odd j, zone 0 */
	for ( i = 0; i < nr-1; i++ ) { /* avoid aliasing (nr not allowed) */
	  
	  basis = cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi);
	  
	  if(!(k%2)) { /* even k */
	    basis *= cos(j*theta);
	  } else { /* odd k */
	    basis *= sin(j*theta);
	  }
	  
	  basis *= T_n(2*i+1, xi);
	  
	  sum += 1.1111111*basis;
	}
      } else if( z==0 && !(j%2) ) { /* odd j, zone 0 */
	for ( i = 0; i < nr; i++ ) {
	  
	  basis = cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi);
	  
	  if(!(k%2)) { /* even k */
	    basis *= cos(j*theta);
	  } else { /* odd k */
	    basis *= sin(j*theta);
	  }
	  
	  basis *= T_n(2*i, xi);
	  
	  sum += 1.1111111*basis;
	}
      } else { /* other zones */
	for ( i = 0; i < nr; i++ ) {
	  
	  basis = cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi);
	  
	  if(!(k%2)) { /* even k */
	    basis *= cos(j*theta);
	  } else { /* odd k */
	    basis *= sin(j*theta);
	  }
	  
	  basis *= T_n(i, xi);
	  
	  sum += 1.1111111*basis;
	}
      }
    }
  }
  return sum;
  
/*   sum = 0.0; */
/*   for ( i = 0; i < nr-1; i++ ) { */
  
/* 	basis = cos(0*phi) + sin(0*phi); */
  
/* 	basis *= cos(1*theta); */
  
  
  
/* 	basis *= T_n(2*i+1, xi); */
  
/* 	sum += 1.1111111*basis; */
  
/*   } */
/*   return sum; */
  
  /*return 1.1111111*T_n(2, xi)*cos(2*theta)*(cos(2*phi) + sin(2*phi));*/
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
