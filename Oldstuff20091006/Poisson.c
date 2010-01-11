/*To compile type: gcc -lm -lgsl -lgslcblas -Wall -pedantic -ansi Poisson.c */
/*To write to file type: ./a.out > phigamma2.txt */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* constants and simple functions */
#define PI 3.141592653589793
#define delta(i, j) (i==j ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) (i%2 ? -1 : 1) /* (-1)^i */
#define cosipiby2(i) (i%2 ? 0 : (i%4==0 ? 1 : -1)) /* \cos(i\pi/2) */

/* function prototypes */

/* the function I ultimately want to create: */
/*void solve_poisson_coeff(tensor_struct_4 *coefficients_ilmn, 
						 tensor_struct_3 *boundary_ilm,
						 tensor_struct_4 *source_ilmn);*/

/* function to evaluate solution at points (r, theta, phi) */
/*double eval_poisson(tensor_struct_4 *coefficients_ilmn, 
					tensor_struct_3 *boundary_ilm, 
					double r, 
					double theta, 
					double phi);*/

double source(double r);

void set_x(gsl_matrix *m);
void set_xinv(gsl_matrix *m);
void set_dbydx(gsl_matrix *m);
void set_xmin1inv(gsl_matrix *m);

void set_A_kernel_even(gsl_matrix *A_even, int L);
void set_A_kernel_odd(gsl_matrix *m, int L);
void set_A_shell(gsl_matrix *m, int L, double alpha, double beta);
void set_A_ext(gsl_matrix *m, int L);

void solve_kernel_even(int L, double alpha, gsl_vector *s_even, gsl_vector *f_even_trunc);
void solve_kernel_odd(int L, double alpha, gsl_vector *s_odd, gsl_vector *f_odd_trunc);
void solve_shell(int L, double alpha, double beta, gsl_vector *s, gsl_vector *f_trunc);
void solve_ext(int L, double alpha, gsl_vector *s, gsl_vector *f_trunc);

void chebfit(double a, double b, gsl_vector *c, double (*func)(double));
double chebeval(double a, double b, gsl_vector *c, double x);

void print_vector(gsl_vector *v);
void print_matrix(gsl_matrix *m);

/* main program */
int main (void)
{
  int N = 9;
  int i, j;
  
  gsl_matrix *some_m = gsl_matrix_alloc (N, N);
  gsl_matrix *A_even = gsl_matrix_alloc (N, N);
  gsl_matrix *A_odd = gsl_matrix_alloc (N, N);
  gsl_matrix *A_shell = gsl_matrix_alloc (N, N);
  gsl_matrix *A_ext = gsl_matrix_alloc (N, N);
  /*gsl_matrix *A_even_trunc = gsl_matrix_alloc (N-1, N-1);*/

  gsl_vector *boundary = gsl_vector_alloc (N);
  gsl_vector *s_kernel = gsl_vector_alloc (2*N-1);
  gsl_vector *s_kernel_even = gsl_vector_alloc (N);
  /*gsl_vector *s_kernel_even_trunc = gsl_vector_alloc (N-1);*/
  gsl_vector *f_trunc = gsl_vector_alloc (N-1);

  double x;

  gsl_vector_set(boundary, 0, 5.0); /* set location of boundary */

  /*set_x(some_m);
  printf("x:\n");
  print_matrix(some_m);

  set_xinv(some_m);
  printf("1/x:\n");
  print_matrix(some_m);

  set_xmin1inv(some_m);
  printf("1/(x-1):\n");
  print_matrix(some_m);*/


  /* Solve system for kernel with a given source */
  /* ------------------------------------------- */
  
  chebfit(-gsl_vector_get(boundary, 0), gsl_vector_get(boundary, 0), s_kernel, source);
  print_vector(s_kernel);
  /* only use even terms */
  for(i=0; i<N; i++)
	gsl_vector_set(s_kernel_even, i, gsl_vector_get(s_kernel, 2*i));
  print_vector(s_kernel_even);

  solve_kernel_even(0, 5, s_kernel_even, f_trunc);

  for(x=0.0; x<=5; x+=0.1)
    printf("%f\t%.18e\n", x, chebeval(-5, 5, f_trunc, x));









  set_A_kernel_odd(A_odd, 1);
  printf("A_kernel_odd:\n");
  print_matrix(A_odd);  
   
 

 
 
  set_A_shell(A_shell, 7, 1.0, 1.0);
  printf("A_shell:\n");
  print_matrix(A_shell);
  
  /* make A_shell matrix banded */
  
  printf("A_shell banded:\n");
  print_matrix(A_shell);
  


  set_A_ext(A_ext, 5);
  printf("A_ext:\n");
  print_matrix(A_ext);

  /* make A_ext matrix banded */
  
  printf("A_ext banded:\n");
  print_matrix(A_ext);



  /* free memory */
  gsl_matrix_free(some_m);
  gsl_matrix_free(A_even);
  gsl_matrix_free(A_odd);
  gsl_matrix_free(A_shell);
  gsl_matrix_free(A_ext);
  
  return 0;
}

/*************************************/
/*            Source                 */
/*************************************/

double source(double r)
{
  double surface = 5.0;
  if(r<=surface)
	return surface - r*r/surface;
  else
	return pow(surface, 5)/pow(r, 4);
}

/*************************************/
/*    Matrix operator functions.     */
/*************************************/

void set_x(gsl_matrix *m)
{
  int n; /*nXn matrix*/
  int i; /*row*/ 
  
  n=(*m).size1;

  /*initialize matrix*/
  gsl_matrix_set_zero(m);
  
  /*set nonzero elements*/
  gsl_matrix_set (m, 0, 1, 0.5);
  gsl_matrix_set (m, 1, 0, 1.0);
  gsl_matrix_set (m, 1, 2, 0.5);
  for(i=2; i<n; i++){
    gsl_matrix_set (m, i, i-1, 0.5);
    if(i+1<n)
      gsl_matrix_set (m, i, i+1, 0.5);
  } 
}

void set_xinv(gsl_matrix *m)
{
  int n;
  int i; /*row*/ 
  int j; /*column*/
  int sign;
  
  n=(*m).size1;

  /*initialize matrix*/
  gsl_matrix_set_zero(m);
 
  /*set nonzero elements*/
  /*first row: (0, 1, 0, -1, 0, 1, 0, -1)*/
  sign=1;
  for(j=1; j<n; j+=2){
    gsl_matrix_set (m, 0, j, sign);
    sign*=-1; /*switch sign*/
  }
  /*all other rows*/
  for(i=1; i<n; i++){
    sign=1;
    for(j=i+1; j<n; j+=2){
      gsl_matrix_set (m, i, j, sign*2);
      sign*=-1;
    }
  }
}


void set_xmin1inv(gsl_matrix *m)
{
  int n; /*nXn matrix*/
  int i; /*row*/
  int j; /*column*/
  
  n=(*m).size1;
  
  /*initialize matrix*/
  gsl_matrix_set_zero(m);
  
  /*set nonzero elements*/
  /*first row*/  
  for(j=1; j<n; j++)
	gsl_matrix_set (m, 0, j, j);
  for(i=1; i<n; i++){
	for(j=i+1; j<n; j++)
	  gsl_matrix_set (m, i, j, 2*(j-i));
  } 
}

/*Returns [0...(n-1)]x[0...(n-1)] matrix for derivative of Chebyshev basis.*/  
void set_dbydx(gsl_matrix *m)
{
  int n;
  int i, j;
  double temp; 
  
  n=(*m).size1;
  
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

/**********************************************************/
/*   make differential operator for even part of kernel   */ 
/*   A =  d^2/dx^2 + (2/x)d/dx - L(L+1)/x^2               */
/**********************************************************/
void set_A_kernel_even(gsl_matrix *A_even, int L)
{
  int i; /*row*/ 
  int j; /*column*/
  int N_even=(*A_even).size1;
  int N_all=2*N_even-1;
  gsl_matrix *dbydx = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *d2bydx2 = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *xinv = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *xinv2 = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *xinvdbydx = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *A = gsl_matrix_calloc (N_all, N_all); /*initialize to 0*/
  
  if(L%2)
	printf("L is not even in the function A_kernel_even.\n");

  /* make differential operator A = d^2/dx^2 + (2/x)d/dx - L(L+1)/x^2 */
  set_dbydx(dbydx);
  set_xinv(xinv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
				 1.0, dbydx, dbydx,
				 0.0, d2bydx2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
				 2.0, xinv, dbydx,
				 0.0, xinvdbydx);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
				 -L*(L+1), xinv, xinv,
				 0.0, xinv2);
  gsl_matrix_add(A, d2bydx2);
  gsl_matrix_add(A, xinvdbydx);
  gsl_matrix_add(A, xinv2);

  /* use only even terms */
  for(i=0; i<N_even; i++)
    for(j=0; j<N_even; j++)
      gsl_matrix_set (A_even, i, j, gsl_matrix_get(A, 2*i, 2*j));
  
  /* free memory */
  gsl_matrix_free(dbydx);
  gsl_matrix_free(d2bydx2);
  gsl_matrix_free(xinv);
  gsl_matrix_free(xinv2);
  gsl_matrix_free(xinvdbydx);
  gsl_matrix_free(A);
}


/**********************************************************/
/*   make differential operator for odd part of kernel    */ 
/*   A =  d^2/dx^2 + (2/x)d/dx - L(L+1)/x^2               */
/**********************************************************/
void set_A_kernel_odd(gsl_matrix *A_odd, int L)
{
  int i; /*row*/ 
  int j; /*column*/
  int N_odd=(*A_odd).size1;
  int N_all=2*N_odd;
  gsl_matrix *dbydx = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *d2bydx2 = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *xinv = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *xinv2 = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *xinvdbydx = gsl_matrix_alloc (N_all, N_all);
  gsl_matrix *A = gsl_matrix_calloc (N_all, N_all); /*initialize to 0*/
  
  if(!(L%2))
	printf("L is not odd in the function A_kernel_odd.\n");

  /* make differential operator A = d^2/dx^2 + (2/x)d/dx - L(L+1)/x^2 */
  set_dbydx(dbydx);
  set_xinv(xinv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
				 1.0, dbydx, dbydx,
				 0.0, d2bydx2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
				 2.0, xinv, dbydx,
				 0.0, xinvdbydx);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
				 -L*(L+1), xinv, xinv,
				 0.0, xinv2);
  gsl_matrix_add(A, d2bydx2);
  gsl_matrix_add(A, xinvdbydx);
  gsl_matrix_add(A, xinv2);

  /* use only even terms */
  for(i=0; i<N_odd; i++)
    for(j=0; j<N_odd; j++)
      gsl_matrix_set (A_odd, i, j, gsl_matrix_get(A, 2*i+1, 2*j+1));
  
  /* free memory */
  gsl_matrix_free(dbydx);
  gsl_matrix_free(d2bydx2);
  gsl_matrix_free(xinv);
  gsl_matrix_free(xinv2);
  gsl_matrix_free(xinvdbydx);
  gsl_matrix_free(A);
}


/**********************************************************/
/*   make differential operator for a shell               */ 
/*   A = (x+b/a)^2(d^2/dx^2) + 2(x+b/a)(d/dx) - L(L+1)    */
/**********************************************************/
void set_A_shell(gsl_matrix *A, int L, double alpha, double beta)
{
  int N=(*A).size1;
  gsl_matrix *id = gsl_matrix_alloc (N, N);
  gsl_matrix *dbydx = gsl_matrix_alloc (N, N);
  gsl_matrix *d2bydx2 = gsl_matrix_alloc (N, N);
  gsl_matrix *x = gsl_matrix_alloc (N, N);
  gsl_matrix *xplusbbya = gsl_matrix_alloc (N, N);
  gsl_matrix *xplusbbya2 = gsl_matrix_alloc (N, N);
  gsl_matrix *temp1 = gsl_matrix_alloc (N, N);
  gsl_matrix *temp2 = gsl_matrix_alloc (N, N);
  gsl_matrix_set_zero(A); /* initialize A to 0 */
  gsl_matrix_set_identity(id); /* initialize A to 0 */
  if(alpha==0.0)
	printf("alpha=0.0 in set_A_shell.  That's a problem.");

  /* make (x+b/a) operator */
  gsl_matrix_set_identity(xplusbbya);
  gsl_matrix_scale(xplusbbya, beta/alpha);
  set_x(x);
  gsl_matrix_add(xplusbbya, x);
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, xplusbbya, xplusbbya,
				 0.0, xplusbbya2);
  
  set_dbydx(dbydx);  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, dbydx, dbydx,
				 0.0, d2bydx2);
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 2.0, xplusbbya, dbydx,
				 0.0, temp1);
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, xplusbbya2, d2bydx2,
				 0.0, temp2);
  
  gsl_matrix_add(A, temp1);
  gsl_matrix_add(A, temp2);
  gsl_matrix_scale(id, -L*(L+1));
  gsl_matrix_add(A, id);
  
  /* free memory */
  gsl_matrix_free(dbydx);
  gsl_matrix_free(d2bydx2);
  gsl_matrix_free(x);
  gsl_matrix_free(xplusbbya);
  gsl_matrix_free(xplusbbya2);
  gsl_matrix_free(temp1);
  gsl_matrix_free(temp2);
}


/**********************************************************/
/*   make differential operator for external domain       */ 
/*   A = d^2/dx^2 - L(L+1)/(x-1)^2                        */
/**********************************************************/
void set_A_ext(gsl_matrix *A, int L)
{
  int N=(*A).size1;
  gsl_matrix *dbydx = gsl_matrix_alloc (N, N);
  gsl_matrix *d2bydx2 = gsl_matrix_alloc (N, N);
  gsl_matrix *xmin1inv = gsl_matrix_alloc (N, N);
  gsl_matrix *xmin1inv2 = gsl_matrix_alloc (N, N);
  
  gsl_matrix_set_zero(A); /* initialize A to 0 */

  /* make differential operator A = d^2/dx^2 - L(L+1)/(x-1)^2 */
  set_dbydx(dbydx);
  set_xmin1inv(xmin1inv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, dbydx, dbydx,
				 0.0, d2bydx2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 -L*(L+1), xmin1inv, xmin1inv,
				 0.0, xmin1inv2);
  gsl_matrix_add(A, d2bydx2);
  gsl_matrix_add(A, xmin1inv2);
  
  /* free memory */
  gsl_matrix_free(dbydx);
  gsl_matrix_free(d2bydx2);
  gsl_matrix_free(xmin1inv);
  gsl_matrix_free(xmin1inv2);
}

/*******************************************************/
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_kernel_even(int L, double alpha, gsl_vector *s_even, gsl_vector *f_even)
{
  int N = (*s_even).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */

  if(L == 0||1)
	p = 1;
  else
	p = 2;
  
  gsl_matrix *A_even = gsl_matrix_alloc (N, N);
  gsl_matrix *A_even_trunc = gsl_matrix_alloc (N-p, N-p);

  gsl_vector *s_even_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_even_trunc = gsl_vector_alloc (N-p);

  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);

  gsl_vector_set_zero(f_even);

  /* make operator matrix */
  set_A_kernel_even(A_even, L);
  
  /* scale s_kernel_even by \alpha^2 */
  gsl_vector_scale(s_even, alpha*alpha);
  
  /* make matrix banded and perform the same operations on source vector */
  /* L_i = (1+\delta_{0i})L_i - L_{i+2}, for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_even, i, j, 
					 (1+delta(i, 0))*gsl_matrix_get(A_even, i, j) - gsl_matrix_get(A_even, i+2, j));
	gsl_vector_set(s_even, i, 
				   (1+delta(i, 0))*gsl_vector_get(s_even, i) - gsl_vector_get(s_even, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_even, i, j, gsl_matrix_get(A_even, i, j) - gsl_matrix_get(A_even, i+2, j));
	gsl_vector_set(s_even, i, 
				   gsl_vector_get(s_even, i) - gsl_vector_get(s_even, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_even, i, j, gsl_matrix_get(A_even, i, j) - gsl_matrix_get(A_even, i+1, j));
	gsl_vector_set(s_even, i, 
				   gsl_vector_get(s_even, i) - gsl_vector_get(s_even, i+1));
  }
  /*printf("A_even banded:\n");
  print_matrix(A_even);
  print_vector(s_even);*/
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_even_trunc, i, j, gsl_matrix_get(A_even, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_even_trunc, i, gsl_vector_get(s_even, i));

  print_matrix(A_even_trunc);
  print_vector(s_even_trunc);
  
  /********************************************************/
  /* TO DO: Use LAPACK to find f instead of what's below, */
  /*        to take advantage of the banded matrix        */
  /********************************************************/

  /* Solve for f_even in A_even_trunc.f_even_trunc = s_even_trunc */
  gsl_linalg_LU_decomp (A_even_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_even_trunc, permute, s_even_trunc, f_even_trunc); 
  print_vector(f_even_trunc);

  for(i=0; i<N-p; i++)
	gsl_vector_set(f_even, i, gsl_vector_get(f_even_trunc, i));

  /* free memory */
  gsl_matrix_free(A_even);
  gsl_matrix_free(A_even_trunc);
  gsl_vector_free(s_even_trunc);
  gsl_vector_free(f_even_trunc);
  gsl_permutation_free(permute);
}

void solve_kernel_odd(int L, double alpha, gsl_vector *s_odd, gsl_vector *f_odd_trunc)
{
 int N = (*s_even).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */

  if(L == 0||1)
	p = 1;
  else
	p = 2;
  
  gsl_matrix *A_even = gsl_matrix_alloc (N, N);
  gsl_matrix *A_even_trunc = gsl_matrix_alloc (N-p, N-p);

  gsl_vector *s_even_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_even_trunc = gsl_vector_alloc (N-p);

  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);

  gsl_vector_set_zero(f_even);

  /* make operator matrix */
  set_A_kernel_even(A_even, L);
  
  /* scale s_kernel_even by \alpha^2 */
  gsl_vector_scale(s_even, alpha*alpha);
  
  /* make matrix banded and perform the same operations on source vector */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+2, j));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+2, j));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+1, j));
  }
  printf("A_kernel_odd banded:\n");
  print_matrix(A_odd);


}


void solve_shell(int L, double alpha, double beta, gsl_vector *s, gsl_vector *f_trunc)
{
 int N = (*s_even).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */

  if(L == 0||1)
	p = 1;
  else
	p = 2;
  
  gsl_matrix *A_even = gsl_matrix_alloc (N, N);
  gsl_matrix *A_even_trunc = gsl_matrix_alloc (N-p, N-p);

  gsl_vector *s_even_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_even_trunc = gsl_vector_alloc (N-p);

  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);

  gsl_vector_set_zero(f_even);

  /* make operator matrix */
  set_A_kernel_even(A_even, L);
  
  /* scale s_kernel_even by \alpha^2 */
  gsl_vector_scale(s_even, alpha*alpha);
  
  /* make matrix banded and perform the same operations on source vector */
  /* L_i = \frac{(1+\delta_{0i})L_i - L_{i+2}}{i+1}, for 0 <= i <= N-3   */
  /* L_i = L_i - L_{i+2},                            for 0 <= i <= N-5   */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_shell, i, j, 
					 ((1+delta(i, 0))*gsl_matrix_get(A_shell, i, j) - gsl_matrix_get(A_shell, i+2, j))/(i+1));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_shell, i, j, gsl_matrix_get(A_shell, i, j) - gsl_matrix_get(A_shell, i+2, j));
  }
}


void solve_ext(int L, double alpha, gsl_vector *s, gsl_vector *f_trunc)
{
 int N = (*s_even).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */

  if(L == 0||1)
	p = 1;
  else
	p = 2;
  
  gsl_matrix *A_even = gsl_matrix_alloc (N, N);
  gsl_matrix *A_even_trunc = gsl_matrix_alloc (N-p, N-p);

  gsl_vector *s_even_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_even_trunc = gsl_vector_alloc (N-p);

  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);

  gsl_vector_set_zero(f_even);

  /* make operator matrix */
  set_A_kernel_even(A_even, L);
  
  /* scale s_kernel_even by \alpha^2 */
  gsl_vector_scale(s_even, alpha*alpha);
  
  /* make matrix banded and perform the same operations on source vector */
  /* L_i = (1+\delta_{0i})L_i - L_{i+2}, for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_ext, i, j, 
					 (1+delta(i, 0))*gsl_matrix_get(A_ext, i, j) - gsl_matrix_get(A_ext, i+2, j));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_ext, i, j, gsl_matrix_get(A_ext, i, j) - gsl_matrix_get(A_ext, i+2, j));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_ext, i, j, gsl_matrix_get(A_ext, i, j) - gsl_matrix_get(A_ext, i+1, j));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_ext, i, j, gsl_matrix_get(A_ext, i, j) - gsl_matrix_get(A_ext, i+2, j));
  }
}


/*******************************************************/
/* Going between functions and Chebyshev coefficients. */
/*******************************************************/


/***********************************************************/
/*Takes a function and finds N Chebyshev coefficients      */
/*upto order N-1 using Chebyshev-Gauss-Lobatto quadrature. */
/***********************************************************/
void chebfit(double a, double b, gsl_vector *c, double (*func)(double))
{
  int N;             /* number of terms in Chebyshev expansion */
  int k;             /* kth coefficient of Chebyshev expansion */
  int j;             /* jth term in summation for coefficient c_k */
  double middle;     /* midpoint of range [a, b] */ 
  double half_width; /* half width of [a, b] */
  double x_j;        /* Gauss-Lobatto quadrature point */
  double sum;        /* sum of terms in summation for c_k */
  double fac;        /* overall factor in front of summation for c_k */
  
  N = (*c).size;
  
  middle = 0.5*(b+a);
  half_width = 0.5*(b-a);
  
  /*Find each of N coefficients c_k.*/
  for(k=0; k<N; k++){
    sum=0;
    /*Evaluate jth term in summation for kth coefficient.*/
    for(j=0; j<N; j++){
      /*Rescale range to be [-1, 1].*/
      x_j = cos(PI*j/(N-1));
      sum += (*func)(half_width*x_j+middle)*cos(PI*k*j/(N-1))/
		(1.0+delta(0, j)+delta(N-1, j));
    }
    fac = 2.0/((N-1)*(1.0+delta(0, k)+delta(N-1, k)));
    gsl_vector_set(c, k, fac*sum);
  }
}

/**************************************************************/
/*Takes N coefficients of Chebyshev polynomials upto order N-1*/
/*and uses Clenshaw's recurrance formula to determine the     */
/*function value at x.                                        */
/**************************************************************/
double chebeval(double a, double b, gsl_vector *c, double x)
{ 
  double t; /*change of variable*/ 
  int N; /*number of polynomials*/
  int k;
  double y_k;
  double y_kplus1; 
  double y_kplus2;  

  N = (*c).size;

  if ((x-a)*(x-b) > 0.0)
    printf("x is not between a and b in chebeval");
  
  /*shift range from [a, b] to [-1, 1]*/
  t = (2.0*x-a-b)/(b-a);
  
  /*Begin Clenshaw's recurrance formula.*/
  y_k = 0.0;
  y_kplus1 = 0.0; 
  y_kplus2 = 0.0;
  for(k=N-1; k>=1; k--){
    y_kplus2 = y_kplus1;
    y_kplus1 = y_k;
    y_k = 2*t*y_kplus1 - y_kplus2 + gsl_vector_get(c, k);
  }

  /*k is now 1*/
  return y_k*t - y_kplus1 + gsl_vector_get(c, 0);
}





/******************************************************/
/*    Functions for printing matrices and vectors.    */
/******************************************************/

void print_vector(gsl_vector *v)
{
  int n;
  int i;

  n = (*v).size;
  
  for (i = 0; i < n; i++){
    printf ("%5g   ", gsl_vector_get (v, i));
  } 
  printf("\n\n");
}


void print_matrix(gsl_matrix *m)
{
  int rows, columns;
  int i, j;

  rows = (*m).size1;
  columns = (*m).size2;
  
  for (i = 0; i < rows; i++){
    for (j = 0; j < columns; j++)
      printf ("%5g   ", gsl_matrix_get (m, i, j));
    printf ("\n");
  } 
  printf("\n");
}
