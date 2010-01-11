/* To compile type: gcc -lm -lfftw3 -lgsl -lgslcblas setgettest.c print.c coefficients.c poisson.h */


/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

int main (void)
{
  int z, i, j, k;
  double xi_i, theta_j, phi_k;
  int nz = 2;
  int nr = 7; /* must be even? */
  int nt = 3;
  int np = 2; /* must be even */
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
	  for ( j = 0; j < nt; j++ ) {
		scalar2d_set(b_zjk, z, j, k, (z+1)*100 + (j+1)*10 + (k+1));
	  }
	}
  }
 
  /* Assign boundary coefficients */
  for ( z = 0; z < nz-1; z++ ) {
	for ( j = 0; j < nt; j++ ) {
	  for ( k = 0; k < npc; k++ ) {
		bound_coeff_set(bcoeff_zjk, z, j, k, REAL, (z+1)*100 + (j+1)*10 + (k+1));
		bound_coeff_set(bcoeff_zjk, z, j, k, IMAG, (z+1)*100 + (j+1)*10 + (k+1));
	  }
	}
  }
  
  /* Assign the field data. */
  for ( z = 0; z < nz; z++ ) {
	for ( i = 0; i < nr; i++ ) {
	  for ( j = 0; j < nt; j++ ) {
		for ( k = 0; k < np; k++ ) {
		  scalar3d_set(c_zijk, z, i, j, k, 
					   (z+1)*1000 + (i+1)*100 + (j+1)*10 + (k+1));
		}
	  }
	}
  }
  
  /* Assign field coefficients */
  for ( z = 0; z < nz; z++ ) {
	for ( i = 0; i < nr; i++ ) {
	  for ( j = 0; j < nt; j++ ) {
		for ( k = 0; k < npc; k++ ) {
		  coeff_set(coeff_zijk, z, i, j, k, REAL, (z+1)*1000 + (i+1)*100 + (j+1)*10 + (k+1));
		  coeff_set(coeff_zijk, z, i, j, k, IMAG, (z+1)*1000 + (i+1)*100 + (j+1)*10 + (k+1));
		}
	  }
	}
  }
  
  /* print start data */
  printf("  Boundary data:\n");  
  print_scalar2d(b_zjk);
  
  printf("  Boundary coefficients:\n");
  print_bound_coeff(bcoeff_zjk);

  printf("  Field data:\n");
  print_scalar3d(c_zijk);

  printf("  Field coefficients:\n");
  print_coeff_2(coeff_zijk);

  return 0;
}
