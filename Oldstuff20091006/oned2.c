/*To compile type: gcc -lm -lgsl -lgslcblas -Wall -pedantic -ansi oned.c */
/*To write to file type: ./a.out > testfile2.txt */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define Ne 10 /*Number of terms in matrices for even functions.*/
#define PI 3.141592653589793
#define delta(i, j) (i==j ? 1 : 0)
#define neg1toi(i) (i%2 ? -1 : 1)
#define cosipiby2(i) (i%2 ? 0 : (i%4==0 ? 1 : -1))

void printvector(gsl_vector *v, int elements);
void printmatrix(gsl_matrix *m, int rows, int columns);
void dchebij(gsl_matrix *m, int n);
void inverse(gsl_matrix *m, int n);
void xmatrix(gsl_matrix *m, int n);
void chebfit(float a, float b, float *c, int N, float (*func)(float));
float source(float x);
float chebeval(float a, float b, float *c, int N, float x);

int
main (void)
{
  int i, j;
  float start=-1.0, end=1.0;
  /*float phisurf=-1.0;*/
  float x;
  int Nall=2*Ne-1;
  float ucoeff[2*Ne-1];
  gsl_matrix *d = gsl_matrix_alloc (Nall, Nall); /*derivative in T_n basis*/
  gsl_matrix *dsq = gsl_matrix_alloc (Nall, Nall); /*2nd derivative operator*/
  gsl_matrix *dsq_even = gsl_matrix_alloc (Ne, Ne);
  gsl_matrix *inv = gsl_matrix_alloc (Nall, Nall); /*1/x operator*/
  gsl_matrix *xop = gsl_matrix_alloc (Nall, Nall); /*x operator*/
  gsl_matrix *L = gsl_matrix_alloc (Nall, Nall); /*linear differential operator*/ 
  gsl_matrix *L_even = gsl_matrix_alloc (Ne, Ne); /*even part of L*/
  gsl_matrix *invD_even = gsl_matrix_alloc (Ne, Ne); /*even part of L*/
  gsl_matrix *T = gsl_matrix_alloc (Ne, Ne); /*T_ki = T_i(x_k)*/
  gsl_matrix *TL = gsl_matrix_alloc (Ne, Ne);
  gsl_vector *s = gsl_vector_alloc (Ne); /*source in T_n basis*/
  gsl_vector *u = gsl_vector_alloc (Ne); /*solution in T_n basis*/
  int k; /*Used in LU decomposition.  I don't know why.*/
  gsl_permutation * p = gsl_permutation_alloc (Ne); /*Used in LU decomposition.*/
  
  
  /*Make linear differential operator L.*/
  /*L = d^2u/dr^2 + (2/r)du/dr*/
  dchebij(d, Nall);
  printmatrix(d, Nall, Nall);
  inverse(inv, Nall);
  printmatrix(inv, Nall, Nall);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, d, d,
		 0.0, dsq);
  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (dsq_even, i, j, gsl_matrix_get(dsq, 2*i, 2*j));
  printf("d^2/dx^2 even part:\n");
  printmatrix(dsq_even, Ne, Ne);

  printf("(1/x)d/dx:\n");
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, inv, d,
		 0.0, L);
  printmatrix(L, Nall, Nall);
  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (invD_even, i, j, gsl_matrix_get(L, 2*i, 2*j));
  printmatrix(invD_even, Ne, Ne);
  


  gsl_matrix_scale(L, 2.0);
  gsl_matrix_add(L, dsq);
  /* gsl_matrix_scale(L, 1/(PI*PI));*/
  printmatrix(L, Nall, Nall);
  
  /*xmatrix(xop, Nall);*/
  /*printmatrix(xop, Nall, Nall);*/

  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (L_even, i, j, gsl_matrix_get(L, 2*i, 2*j));
  printmatrix(L_even, Ne, Ne);
  
  /*Make colocation matrix T_ki = T_i(x_k) where x_k = sin(k*PI/(2*ORDER)).*/
  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (T, i, j, neg1toi(j)*cos(i*j*PI/(Ne-1)));
  printmatrix(T, Ne, Ne);
  
  /*Make operator T.L*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, T, L_even,
		 0.0, TL);
  printmatrix(TL, Ne, Ne);
  
  /*Make source vector.*/
  for(i=0; i<Ne; i++){
    gsl_vector_set(s, i, source(sin(i*PI/(2*(Ne-1)))));
    printf("%f\n", sin(i*PI/(2*(Ne-1))));
  }
  printvector(s, Ne);
  
  /*Set boundary conditions in TL and s.*/
  for(i=0; i<Ne; i++){
    gsl_matrix_set (TL, 0, i, neg1toi(i));
  }
  gsl_vector_set(s, 0, -2);
  printmatrix(TL, Ne, Ne);
  printvector(s, Ne);
  
  /*Solve for u~ in T.L.u~ = s~.*/
  gsl_linalg_LU_decomp (TL, p, &k); 
  gsl_linalg_LU_solve (TL, p, s, u); 
  printvector(u, Ne);
  
  ucoeff[0] = gsl_vector_get(u, 0);
  for(i=1; i<Ne; i++){
    ucoeff[2*i-1] = 0;
    ucoeff[2*i] = gsl_vector_get(u, i);
  }
  
  for(x=0.0; x<=end; x+=0.01)
    printf("%f\t%f\n", x, chebeval(start, end, &ucoeff[0], Nall-1, x));
  
  
  gsl_matrix_free(d);
  gsl_matrix_free(dsq);
  gsl_matrix_free(L);
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
/*RECHECK THIS MATRIX.*/
void xmatrix(gsl_matrix *m, int n)
{
  int i; /*row*/ 
  
  /*initialize matrix*/
  gsl_matrix_set_zero(m);
  
  /*set nonzero elements*/
  gsl_matrix_set (m, 0, 1, 1.0);
  for(i=1; i<n; i++){
    gsl_matrix_set (m, i, i-1, 0.5);
    if(i+1<n)
      gsl_matrix_set (m, i, i+1, 0.5);
  } 
}

void inverse(gsl_matrix *m, int n)
{
  int i; /*row*/ 
  int j; /*column*/
  int sign;
  
  /*initialize matrix*/
  gsl_matrix_set_zero(m);
 
  /*set nonzero elements*/
  /*first row: (0, 1, 0, -1, 0, 1, 0, -1)*/
  sign=1;
  for(j=1; j<n; j+=2){
    gsl_matrix_set (m, 0, j, sign);
    sign*=-1; /*switch sign*/
  }
  for(i=1; i<n; i++){
    sign=1;
    for(j=i+1; j<n; j+=2){
      gsl_matrix_set (m, i, j, sign*2);
      sign*=-1;
    }
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

float source(float xi)
{
  if(xi<0.01)
    return PI*PI*(1-pow(xi*PI, 2)/6.0+pow(xi*PI, 4)/120.0);
  else
    return PI*PI*sin(PI*xi)/(PI*xi);
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
