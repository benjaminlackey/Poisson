/*To compile type: gcc -lm -lgsl -lgslcblas -Wall -pedantic -ansi pseudospectral.c */
/*To write to file type: ./a.out > testfile2.txt */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define ORDER 6 /*Highest Chebyshev polynomial. Should never be 0.*/
#define PI 3.141592653589793
#define delta(i, j) (i==j ? 1 : 0)
#define neg1toi(i) (i%2 ? -1 : 1)

void printvector(gsl_vector *v, int elements);
void printmatrix(gsl_matrix *m, int rows, int columns);
void makematrix(gsl_matrix *m, int n);
void dchebij(gsl_matrix *m, int n);
void chebfit(float a, float b, float *c, int N, float (*func)(float));
float source(float x);
float chebeval(float a, float b, float *c, int N, float x);

int
main (void)
{

  int n=ORDER+1; /*number of rows/columns in matrices*/
  int i, j;
  float start=-1.0, end=1.0;
  float x;
  float ucoeff[ORDER+1];
  gsl_matrix *d = gsl_matrix_alloc (n, n); /*derivative in T_n basis*/
  gsl_matrix *dsq = gsl_matrix_alloc (n, n); /*2nd derivative operator*/
  gsl_matrix *L = gsl_matrix_alloc (n, n); /*linear differential operator*/ 
  gsl_matrix *T = gsl_matrix_alloc (n, n); /*T_ki = T_i(x_k)*/
  gsl_matrix *TL = gsl_matrix_alloc (n, n);
  gsl_vector *s = gsl_vector_alloc (n); /*source in T_n basis*/
  gsl_vector *u = gsl_vector_alloc (n); /*solution in T_n basis*/
  int k; /*Used in LU decomposition.  I don't know why.*/
  gsl_permutation * p = gsl_permutation_alloc (n); /*Used in LU decomposition.*/


  /*Make linear differential operator L.*/
  /*L = d^2u/dx^2 - 4du/dx + 4u*/
  dchebij(d, n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, d, d,
		 0.0, dsq);
  gsl_matrix_set_identity(L);
  gsl_matrix_sub(L, d);
  gsl_matrix_scale(L, 4.0);
  gsl_matrix_add(L, dsq);
  /*printmatrix(L, n, n);*/
  
  /*Make colocation matrix T_ki = T_i(x_k) where x_k = cos(k*PI/ORDER).*/
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      gsl_matrix_set (T, i, j, cos(i*j*PI/ORDER));
  /*printmatrix(T, n, n);*/
  
  /*Make operator T.L*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, T, L,
		 0.0, TL);
  /*printmatrix(TL, n, n);*/
  
  /*Make source vector.*/
  for(i=0; i<n; i++)
    gsl_vector_set(s, i, source(cos(i*PI/ORDER)));
  /*printvector(s, n);*/
  
  /*Set boundary conditions in TL and s.*/
  for(i=0; i<n; i++){
    gsl_matrix_set (TL, 0, i, 1);
    gsl_matrix_set (TL, n-1, i, neg1toi(i));
  }
  gsl_vector_set(s, 0, 0);
  gsl_vector_set(s, n-1, 0);
  /*printmatrix(TL, n, n);*/

  /*Solve for u~ in T.L.u~ = s~.*/
  gsl_linalg_LU_decomp (TL, p, &k); 
  gsl_linalg_LU_solve (TL, p, s, u); 
  /*printvector(u, n);*/

    
  for(i=0; i<n; i++){
    ucoeff[i] = gsl_vector_get(u, i);
  }
  
  for(x=start; x<=end; x+=0.1)
    printf("%f\t%f\n", x, chebeval(start, end, &ucoeff[0], ORDER, x));
  
  gsl_matrix_free (d);
  gsl_matrix_free (dsq);
  gsl_matrix_free (L);
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


void printmatrix(gsl_matrix *m, int rows, int columns)
{
  int i, j;
  for (i = 0; i < rows; i++){
    for (j = 0; j < columns; j++)
      printf ("%4g  ", gsl_matrix_get (m, i, j));
    printf ("\n");
  } 
  printf("\n");
}

void printvector(gsl_vector *v, int elements)
{
  int i;
  for (i = 0; i < elements; i++){
    printf ("%4g  ", gsl_vector_get (v, i));
  } 
  printf("\n\n");
}

float source(float x)
{
  /*return pow(x, 4)*exp(-x*x);*/
  return exp(x)-4*exp(1)/(1+exp(2));
  /*return 8*x*x*x*x + 4*x*x*x - 6*x*x - 2*x + 2;*/
}

void chebfit(float a, 
	     float b, 
	     float *c, 
	     int N, 
	     float (*func)(float))
     /*Takes a function and finds the Chebyshev coefficients    *
      *upto order N using Chebyshev-Gauss-Lobatto quadrature.   *
      *There are N+1 coefficients so c should have N+1 elements.*/
{
  int k;            /*kth coefficient of Chebyshev expansion*/
  int j;            /*jth term in summation for coefficient c[k]*/
  float middle;     /*midpoint of range [a, b]*/ 
  float half_width; /*half width of [a, b]*/
  float x_j;        /*Gauss-Lobatto quadrature point*/
  double sum;       /*sum of terms in summation for c[k]*/
  float fac;        /*overall factor in front of summation for c[k]*/
  
  middle = 0.5*(b+a);
  half_width = 0.5*(b-a);
  
  /*Find each of N+1 coefficients c[k] in range [0, N].*/
  for(k=0; k<=N; k++){
    sum=0;
    /*Evaluate jth term in summation for kth coefficient.*/
    for(j=0; j<=N; j++){
      /*Rescale range to be [-1, 1].*/
      x_j = cos(PI*j/N);
      sum += (*func)(half_width*x_j+middle)*cos(PI*k*j/N)/(1.0+delta(0, j)+delta(N, j));
      /*printf("%f\n", sum);*/
    }
    fac = 2.0/(N*(1.0+delta(0, k)+delta(N, k)));
    /*printf("fac is %f\n", fac);*/
    c[k] = fac*sum;
    /*printf("c[%d] is %f\n", k, c[k]);*/
  }
}

float chebeval(float a, 
	       float b, 
	       float *c, 
	       int N, 
	       float x)
     /*Takes coefficients of Chebyshev polynomials upto order N*
      *and uses Clenshaw's recurrance formula to determine the *
      *function value at x.                                    */
{ 
  float t; /*change of variable*/ 
  int k;
  float y_k;
  float y_kplus1; 
  float y_kplus2;  
  
  if ((x-a)*(x-b) > 0.0)
    printf("x is not between a and b in chebeval");
  
  /*shift range from [a, b] to [-1, 1]*/
  t = (2.0*x-a-b)/(b-a);
  
  /*Begin Clenshaw's recurrance formula.*/
  y_k = 0.0;
  y_kplus1 = 0.0; 
  y_kplus2 = 0.0;
  for(k=N; k>=1; k--){
    y_kplus2 = y_kplus1;
    y_kplus1 = y_k;
    y_k = 2*t*y_kplus1 - y_kplus2 + c[k];
  }

  /*k is now 1*/
  return y_k*t - y_kplus1 + c[0];
}
