/*To compile type: gcc -lm -lgsl -lgslcblas -Wall -pedantic -ansi test3.c */
/*To write to file type: ./a.out > testfile2.txt */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define ORDER 4 /*Highest Chebyshev polynomial*/
#define PI 3.141592653589793
#define delta(i, j) (i==j ? 1 : 0)
#define min1toint(i) (i%2 ? -1 : 1)

void printmatrix(gsl_matrix *m, int n);
void makematrix(gsl_matrix *m, int n);
void dchebij(gsl_matrix *m, int n);
void chebfit(float a, float b, float *c, int n, float (*func)(float));
float source(float x);
float chebeval(float a, float b, float *c, int m, float x);

int
main (void)
{

  int n=ORDER+1; /*number of rows/columns in matrices*/
  int i, j;
  float start=-1.0, end=1.0;
  float y;
  float coeff[ORDER+1];
  float ucoeff[ORDER+1];
  gsl_matrix *d = gsl_matrix_alloc (n, n);
  gsl_matrix *dsq = gsl_matrix_alloc (n, n);
  gsl_matrix *op = gsl_matrix_alloc (n, n);
  gsl_vector *s = gsl_vector_alloc (n);
  
  dchebij(d, n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, d, d,
		 0.0, dsq);
  gsl_matrix_set_identity(op);
  
  gsl_matrix_sub(op, d);
  gsl_matrix_scale(op, 4.0);
  gsl_matrix_add(op, dsq);
  
  for(j=0; j<n; j++){
    gsl_matrix_set (op, n-1, j, 1.0);
    gsl_matrix_set (op, n-2, j, min1toint(j));
  }
  
  /*printmatrix(op, n);*/
  
  chebfit(start, end, &coeff[0], n, source);
  /*set boundary conditions for source*/
  coeff[n-1]=0;
  coeff[n-2]=0;
  
  for(i=0; i<n-2; i++)
    gsl_vector_set(s, i, coeff[i]);
  
  for(i=0; i<n; i++) 
    printf("%f\n", gsl_vector_get(s, i)); 
  
  
  gsl_vector *x = gsl_vector_alloc (n);
  
  int k;
  
  gsl_permutation * p = gsl_permutation_alloc (n);
  
  gsl_linalg_LU_decomp (op, p, &k);
  
  gsl_linalg_LU_solve (op, p, s, x);
  
  printf ("x = \n");
  gsl_vector_fprintf (stdout, x, "%g");
  
  for(i=0; i<n; i++)
    ucoeff[i] = gsl_vector_get(x, i);
  
  for(y=start; y<=end; y+=0.1) 
    printf("%f\t%f\n", y, chebeval(start, end, &ucoeff[0], n-1, y)); 
  
  
  
  
  gsl_matrix_free (d);
  gsl_matrix_free (dsq);
  gsl_matrix_free (op);
  return 0;
}

/*Returns [0...(n-1)]x[0...(n-1)] matrix for derivative of Chebyshev basis.*/  
void dchebij(gsl_matrix *m, int n)
{
  int i, j;
  float temp; 
  
  /*initialize matrix*/
  gsl_matrix_set_zero(m);
 
  /*evaluate nonzero elements*/
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 2, 4.0);
  for(i=0; i<n; i++)
    for(j=((i>3) ? i : 3); j<n; j++){
      temp = 2*j*delta(i, j-1) + (j*gsl_matrix_get(m, i, j-2))/(j-2);
      gsl_matrix_set (m, i, j, temp);
    }
}


void printmatrix(gsl_matrix *m, int n)
{
  int i, j;
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++)
      printf ("%4g  ", gsl_matrix_get (m, i, j));
    printf ("\n");
  } 
}

float source(float x)
{
  /*return pow(x, 4)*exp(-x*x);*/
  return exp(x)-4*exp(1)/(1+exp(2));
}

void chebfit(float a, float b, float *c, int n, float (*func)(float))
{
  /*printf("Hello world.");*/
  int k, j;
  float fac, bpa, bma, *f;
  
  f = malloc(sizeof(float)*n);
  
  bma = 0.5*(b-a); /*Width of right side of range.*/
  bpa = 0.5*(b+a); /*Midpoint of range.*/
  
  for(k=0; k<n; k++){
    float zero_scaled = cos(PI*(k+0.5)/n); /*Find kth zero of T_n on [-1, 1]*/
    float zero = zero_scaled*bma + bpa; 
    f[k] = (*func)(zero);
    /*f[k] = exp(-zero)+zero;*/ /*Function value at each zero.*/
  }
  
  fac = 2.0/n; /*Factor in front of summation term.*/
  
  /*Evaluate coefficient c[j] for each T_j.*/
  for(j=0; j<n; j++){
    double sum = 0.0;
    for(k=0; k<n; k++)
      sum += f[k]*cos(PI*j*(k+0.5)/n);
    c[j] = fac*sum;
    /*printf("%f", c[j]);*/
  }
  free(f);
}	

float chebeval(float a, float b, float *c, int m, float x)
{
  float d=0.0, dd=0.0, sv, y, y2;
  int j;
  
  if ((x-a)*(x-b) > 0.0)
    printf("x is not between a and b");
  
  y = (2.0*x-a-b)/(b-a);
  y2 = 2.0*y;
 
  for(j=m-1; j>=1; j--){
    sv = d;
    d = y2*d-dd+c[j];
    dd = sv;
  }

  return y*d-dd+0.5*c[0];
}
