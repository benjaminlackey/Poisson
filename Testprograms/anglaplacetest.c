/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/poisson.h anglaplacetest.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* fftw header */
#include <fftw3.h>

/* own header */
#include "poisson.h"

/*double boundary(int z, int nt, int np, double theta, double phi);*/
double field(int z, int nr, int nt, int np, double xi, double theta, double phi);
double fieldphysical(int z, double r, double theta, double phi);
double T_n(int n, double x);

int main (void)
{
  int z, imag, i, l, m;
  int nz = 1;
  int nr = 1; /* must be odd? */
  int nt = 17;
  int np = 32; /* must be even */
  int npc;
  gsl_matrix **fouriertoylm;
  gsl_matrix **ylmtofourier;
  ylm_coeff *c_ylm;
  ylm_coeff *lapc_ylm;
  coeff *c_fourier;
  coeff *lapc_fourier;
  
  npc = np/2 + 1;
  
  /*>>>>>>>>>>>>>>>> make Fourier to Y_lm matrices <<<<<<<<<<<<<<<*/
  
  fouriertoylm = fouriertoylm_matrix_alloc(nt);
  fouriertoylm_matrix_set(fouriertoylm);

/*   fouriertoylm = malloc(sizeof(gsl_matrix *)*npc); /\* the size of a gsl_matrix pointer times the number npc *\/ */
/*   /\* set even matrices *\/ */
/*   for(m=0; m<nt; m += 2) { */
/*     fouriertoylm[m] = gsl_matrix_alloc(nt-m, nt); */
/*     set_fouriertoylm_meven(fouriertoylm[m], m, nt); */
/*   } */
/*   /\* set odd matrices *\/ */
/*   for(m=1; m<nt; m += 2) { */
/*     fouriertoylm[m] = gsl_matrix_alloc(nt-m, nt); */
/*     set_fouriertoylm_modd(fouriertoylm[m], m, nt); */
/*   } */
  
/*   for(m=0; m<nt; m++) { */
/*     printf("m = %d:\n", m); */
/*     print_matrix(fouriertoylm[m]); */
/*   } */
  
  /*>>>>>>>>>>>>>>>> make Y_lm to Fourier matrices <<<<<<<<<<<<<<<*/
  
  ylmtofourier = ylmtofourier_matrix_alloc(nt);
  ylmtofourier_matrix_set(ylmtofourier);

/*   ylmtofourier = malloc(sizeof(gsl_matrix *)*np); */
/*   /\* set even matrices *\/ */
/*   for(m=0; m<nt; m += 2) { */
/*     ylmtofourier[m] = gsl_matrix_alloc(nt, nt-m); */
/*     set_ylmtofourier_meven(ylmtofourier[m], m, nt); */
/*   } */
/*   /\* set odd matrices *\/ */
/*   for(m=1; m<nt; m += 2) { */
/*     ylmtofourier[m] = gsl_matrix_alloc(nt, nt-m); */
/*     set_ylmtofourier_modd(ylmtofourier[m], m, nt); */
/*   } */
  
/*   for(m=0; m<nt; m++) { */
/*     printf("m = %d:\n", m); */
/*     print_matrix(ylmtofourier[m]); */
/*   } */
  
  c_fourier = coeff_alloc(nz, nr, nt, np); 
  lapc_fourier = coeff_alloc(nz, nr, nt, np);
  c_ylm = ylm_coeff_alloc(nz, nr, nt, npc);
  lapc_ylm = ylm_coeff_alloc(nz, nr, nt, npc);

  /* Assign the spherical harmonic coefficients. */
  for ( z = 0; z < nz; z++ ) {
    for ( imag=0; imag<=1; imag++ ) {
      for ( i = 0; i < nr; i++ ) {
	for ( m = imag; m < nt-imag; m++ ) {
	  for ( l = m; l < nt-m%2; l++ ) {
	    ylm_coeff_set(c_ylm, z, i, l, m, imag, 1.0);
	  }
	}
      }
    }
  }
  print_ylm_coeff(c_ylm);
  
  /* go to fourier basis */
  transform_ylmtofourier(c_ylm, c_fourier, ylmtofourier);
  print_coeff(c_fourier);
 
/*   coeff_set(c_fourier, 0, 0, 0, 4, REAL, 3.0/8.0); */
/*   coeff_set(c_fourier, 0, 0, 2, 4, REAL, -1.0/2.0); */
/*   coeff_set(c_fourier, 0, 0, 4, 4, REAL, 1.0/8.0); */
/*   print_coeff(c_fourier); */

  /* take angular laplacian of c_fourier and return lapc_fourier */
  laplace_ang(lapc_fourier, c_fourier);
  print_coeff(lapc_fourier);
  
  /* return to spherical harmonics */
  transform_fouriertoylm(lapc_fourier, lapc_ylm, fouriertoylm);
  print_ylm_coeff(lapc_ylm);
  

  fouriertoylm_matrix_free(fouriertoylm);
  ylmtofourier_matrix_free(ylmtofourier);

/*  /\* free matrices as well as the array of pointers to the matrices *\/ */
/*   for(m=0; m<nt; m++) { */
/* 	gsl_matrix_free(fouriertoylm[m]); */
/*   } */
/*   free(fouriertoylm); */

/*   /\* free matrices as well as the array of pointers to the matrices *\/ */
/*   for(m=0; m<nt; m++) { */
/* 	gsl_matrix_free(ylmtofourier[m]); */
/*   } */
/*   free(ylmtofourier); */

  return 0;
}
