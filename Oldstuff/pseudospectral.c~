/*To compile type: gcc -lm -lgsl -lgslcblas -Wall -pedantic -ansi galerkin.c */
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
  float coeff[ORDER+1];
  float ucoeff[ORDER+1];
  gsl_matrix *d = gsl_matrix_alloc (n, n); /*derivative in T_n basis*/
  gsl_matrix *dsq = gsl_matrix_alloc (n, n); /*2nd derivative operator*/
  gsl_matrix *op = gsl_matrix_alloc (n, n); /*linear differential operator*/
  gsl_vector *s = gsl_vector_alloc (n); /*source in T_n basis*/
  gsl_vector *s_in_phi = gsl_vector_alloc (n-2); /*source in phi_n basis*/ 
  gsl_matrix *phi = gsl_matrix_calloc (n, n-2); /*transform from T_n to phi_n basis*/
  gsl_matrix *phit_l_phi = gsl_matrix_alloc (n-2, n-2);
  gsl_matrix *temp = gsl_matrix_alloc (n, n-2); 
  gsl_vector *u_in_phi = gsl_vector_alloc (n-2); /*solution in phi_n basis*/
  gsl_vector *u = gsl_vector_alloc (n); /*solution in T_n basis*/
  int k; /*Used in LU decomposition.  I don't know why.*/
  gsl_permutation * p = gsl_permutation_alloc (n-2); /*Used in LU decomposition.*/


  /*make operator op*/
  dchebij(d, n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, d, d,
		 0.0, dsq);
  gsl_matrix_set_identity(op);
  
  gsl_matrix_sub(op, d);
  gsl_matrix_scale(op, 4.0);
  gsl_matrix_add(op, dsq);
  
  /*make transformation matrix phi to go to Galerkin basis*/
  for(j=0; j<n-2; j++){
    /*even elements in 0th row are -1 [-1 0 -1 0 -1...]*/
    if(!(j%2))
       gsl_matrix_set (phi, 0, j, -1.0);
    /*odd elements in 1st row are -1 [0 -1 0 -1 0...]*/
    if(j%2)
      gsl_matrix_set (phi, 1, j, -1.0);
    /*displaced identity matrix*/
    gsl_matrix_set (phi, j+2, j, 1.0);
  }
  /*printmatrix(phi, n, n-2);*/
  
  /*make operator phi^T.L.phi*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, op, phi,
		 0.0, temp);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 
		 1.0, phi, temp,
		 0.0, phit_l_phi);
  /*printmatrix(phit_l_phi, n-2, n-2);*/
  
  /*make source basis and transform to Galerkin basis*/
  
  chebfit(start, end, &coeff[0], ORDER, source); 
  /*for(i=0; i<n; i++) 
    printf("%f\n", coeff[i]);*/
  
  for(i=0; i<n-2; i++)
    gsl_vector_set(s, i, coeff[i]);
  
  /*for(i=0; i<n; i++) 
    printf("%f\n", gsl_vector_get(s, i)); */
  
  gsl_blas_dgemv(CblasTrans, 
		 1.0, phi, s,
		 0.0, s_in_phi);
  /*printvector(s_in_phi, n-2);*/
  
  /*Solve for u~~ in phi^T.L.phi.u~~ = s~~.*/
  gsl_linalg_LU_decomp (phit_l_phi, p, &k);
  gsl_linalg_LU_solve (phit_l_phi, p, s_in_phi, u_in_phi);
  
  /*Transform u back to T_n basis using u~ = phi^T.u~~.*/
  gsl_blas_dgemv(CblasNoTrans, 
		 1.0, phi, u_in_phi,
		 0.0, u);
    
  for(i=0; i<n; i++){
    ucoeff[i] = gsl_vector_get(u, i);
  }
  
  for(x=start; x<=end; x+=0.01) 
    printf("%f\t%f\n", x, chebeval(start, end, &ucoeff[0], ORDER, x)); 
  
  
  
  
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
