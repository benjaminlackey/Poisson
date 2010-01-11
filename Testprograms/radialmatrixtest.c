/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/remainder.c /Users/lackey/Research/Poisson/poisson.h radialmatrixtest.c */

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

int main (void)
{
  int z, L;
  int nz = 3;
  int nr = 9;
  int nt = 5;
  gsl_vector *alpha_v;
  gsl_vector *beta_v;
  gsl_matrix ***radial_matrices;
  
  alpha_v = gsl_vector_alloc(nz);
  beta_v = gsl_vector_alloc(nz);
  radial_matrices = radial_matrix_alloc(nz, nr, nt);

  gsl_vector_set(alpha_v, 0, 5.0);
  gsl_vector_set(alpha_v, 1, 5.0);
  gsl_vector_set(alpha_v, 2, -2.0);
  gsl_vector_set(beta_v, 0, 0.0);
  gsl_vector_set(beta_v, 1, 10.0);
  gsl_vector_set(beta_v, 2, 0.0);
  
  radial_matrix_set(nz, nt, alpha_v, beta_v, radial_matrices);

  for(z = 0; z < nz; z++) {
    for(L = 0; L < nt; L++) {
      printf("z = %d, L = %d\n", z, L); 
      print_matrix(radial_matrices[z][L]);
    }
  }

  gsl_vector_free(alpha_v);
  gsl_vector_free(beta_v);
  radial_matrix_free(nz, nt, radial_matrices);

  return 0;
}
