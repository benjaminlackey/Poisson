/*To compile type: gcc -lm -lgsl -lgslcblas -Wall -pedantic test2.c */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>

void makematrix(gsl_matrix *m, int n);

int
main (void)
{
  int i, j;
  int n=5;
  
  gsl_matrix *m = gsl_matrix_alloc (n, n);
    
  makematrix(m, n);
  
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      printf ("m(%d,%d) = %g\n", 
	      i, j, gsl_matrix_get (m, i, j));
  
  gsl_matrix_free (m);
  
  return 0;
}

void makematrix(gsl_matrix *m, int n)
{
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j <n; j++)
      gsl_matrix_set (m, i, j, 0.23 + 100*i + j);
}


