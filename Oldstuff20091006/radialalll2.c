/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W radialalll2.c */
/* To write to file type: ./a.out > phigamma2.txt */

/*****************************************************************************/
/* TO DO: Make chebfit_even chebfit_odd chebeval_even chebeval_odd functions */
/*****************************************************************************/

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
double source(double r);
double source_ext(double u);

void set_x(gsl_matrix *m);
void set_xinv(gsl_matrix *m);
void set_dbydx(gsl_matrix *m);
void set_xmin1inv(gsl_matrix *m);

void set_A_kernel_even(gsl_matrix *A_even, int L);
void set_A_kernel_odd(gsl_matrix *m, int L);
void set_A_shell(gsl_matrix *m, int L, double alpha, double beta);
void set_A_ext(gsl_matrix *m, int L);

void solve_kernel_even(int L, double alpha, gsl_vector *s_even, gsl_vector *f_even);
void solve_kernel_odd(int L, double alpha, gsl_vector *s_odd, gsl_vector *f_odd);
void solve_shell(int L, double alpha, double beta, gsl_vector *s_old, gsl_vector *f);
void solve_ext(int L, double alpha, gsl_vector *s_old, gsl_vector *f);

void solve_continuity(int L, 
					  gsl_vector *boundary, 
					  gsl_matrix *particular_m, 
					  gsl_vector *homo_coeff);

void chebfit(double a, double b, gsl_vector *c, double (*func)(double));
double chebeval(double a, double b, gsl_vector *c, double x);

void print_vector(gsl_vector *v);
void print_matrix(gsl_matrix *m);

void solve_radial(int L, 
				  int N,
				  double (*source)(double), 
				  double (*source_ext)(double), 
				  gsl_vector *boundary, 
				  gsl_matrix *particular_m, 
				  gsl_vector *homo_coeff); 

double radialeval(int L, 
				  gsl_vector *boundary, 
				  gsl_matrix *particular_m, 
				  gsl_vector *homo_coeff, 
				  double r);

/* main program */
int main (void)
{
  int N = 5;
  int L = 0;
  
  int nz = 2; /* number of zones (must be >= 2) */
  gsl_vector *homo_coeff = gsl_vector_alloc (2*nz-2);
  gsl_vector *boundary = gsl_vector_alloc (nz-1); /* position of each boundary */
  /*                |Chebyshev coefficients of solution in zone 0   | */
  /*  		    |Chebyshev coefficients of solution in zone 1   | */
  /* particular_m = |                       ...                     | */
  /*                |Chebyshev coefficients of solution in zone nz-1| */
  gsl_matrix *particular_m = gsl_matrix_alloc (nz, N);

  gsl_vector_set(boundary, 0, 1.0); /* set location of boundary */

  double r;

  solve_radial(L, N, source, source_ext, boundary, particular_m, homo_coeff); 
 
  /*print_vector(boundary);
  print_matrix(particular_m);
  print_vector(homo_coeff);*/
  
  for(r = 0.0; r <=10.0; r += 0.1)
    printf("%f\t%.18e\n", r, radialeval(L, boundary, particular_m, homo_coeff, r));
  
  /* free memory */
  gsl_vector_free(homo_coeff);
  gsl_vector_free(boundary);
  gsl_matrix_free(particular_m);
  
  return 0;
}

/*************************************/
/*            Source                 */
/*************************************/

/* double source(double r) */
/* { */
/*   double surface = 5.0; */
/*   if(r<=surface) */
/* 	return surface - r*r/surface; */
/*   else */
/* 	return pow(surface, 5)/pow(r, 4); */
/* } */
double source(double r)
{
  double surface = 1.0;
  return surface - r*r/surface;
}

double source_ext(double u)
{
  double surface = 1.0;
  return pow(surface, 5)*pow(u, 4);
}

 
/*************************************/
/*    Matrix operator functions.     */
/*************************************/

/****************************************************/
/* Matrix operator for x acting on Chebyshev basis. */
/* Returns [0...(N-1)]x[0...(N-1)] matrix,          */
/* where N-1 is highest order polynomial.           */
/****************************************************/
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


/******************************************************/
/* Matrix operator for 1/x acting on Chebyshev basis. */
/* Returns [0...(N-1)]x[0...(N-1)] matrix,            */
/* where N-1 is highest order polynomial.             */
/******************************************************/
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


/**********************************************************/
/* Matrix operator for 1/(x-1) acting on Chebyshev basis. */
/* Returns [0...(N-1)]x[0...(N-1)] matrix,                */
/* where N-1 is highest order polynomial.                 */
/**********************************************************/
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


/*******************************************************/
/* Matrix operator for d/dx acting on Chebyshev basis. */
/* Returns [0...(N-1)]x[0...(N-1)] matrix,             */
/* where N-1 is highest order polynomial.              */
/*******************************************************/
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
  gsl_matrix_set_identity(id); /* initialize A to identity */
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
/* Solve even part of kernel:                          */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_kernel_even(int L, double alpha, gsl_vector *s_even, gsl_vector *f_even)
{
  int N = (*s_even).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */
  
  if(L%2)
	printf("L is not even in solve_kernel_even.");
  else if(L == 0)
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
  /* and last p rows of source vector                      */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_even_trunc, i, j, gsl_matrix_get(A_even, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_even_trunc, i, gsl_vector_get(s_even, i));
  /*printf("A_even banded truncated:\n");
  print_matrix(A_even_trunc);
  print_vector(s_even_trunc);*/
  
  /********************************************************/
  /* TO DO: Use LAPACK to find f instead of what's below, */
  /*        to take advantage of the banded matrix        */
  /********************************************************/

  /* Solve for f_even in A_even_trunc.f_even_trunc = s_even_trunc */
  gsl_linalg_LU_decomp (A_even_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_even_trunc, permute, s_even_trunc, f_even_trunc); 
  /*print_vector(f_even_trunc);*/

  /* Set solution vector f_even = (0, ..., f_even_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_even, i+p, gsl_vector_get(f_even_trunc, i));

  /* free memory */
  gsl_matrix_free(A_even);
  gsl_matrix_free(A_even_trunc);
  gsl_vector_free(s_even_trunc);
  gsl_vector_free(f_even_trunc);
  gsl_permutation_free(permute);
}


/*******************************************************/
/* Solve odd part of kernel:                           */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_kernel_odd(int L, double alpha, gsl_vector *s_odd, gsl_vector *f_odd)
{
  int N = (*s_odd).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */
  
  if(!(L%2))
	printf("L is not odd in solve_kernel_odd.");
  else if(L == 1)
	p = 1;
  else
	p = 2;
  
  gsl_matrix *A_odd = gsl_matrix_alloc (N, N);
  gsl_matrix *A_odd_trunc = gsl_matrix_alloc (N-p, N-p);
  
  gsl_vector *s_odd_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_odd_trunc = gsl_vector_alloc (N-p);

  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);

  gsl_vector_set_zero(f_odd);

  /* make operator matrix */
  set_A_kernel_even(A_odd, L);
  
  /* scale s_kernel_even by \alpha^2 */
  gsl_vector_scale(s_odd, alpha*alpha);
  
  /* make matrix banded and perform the same operations on source vector */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+2, j));
	gsl_vector_set(s_odd, i, 
				   gsl_vector_get(s_odd, i) - gsl_vector_get(s_odd, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+2, j));
	gsl_vector_set(s_odd, i, 
				   gsl_vector_get(s_odd, i) - gsl_vector_get(s_odd, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j, gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+1, j));
	gsl_vector_set(s_odd, i, 
				   gsl_vector_get(s_odd, i) - gsl_vector_get(s_odd, i+1));
  }
  /*printf("A_kernel_odd banded:\n");
	print_matrix(A_odd);	
	print_vector(s_odd);*/
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_odd_trunc, i, j, gsl_matrix_get(A_odd, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_odd_trunc, i, gsl_vector_get(s_odd, i));
  
  /*print_matrix(A_odd_trunc);
	print_vector(s_odd_trunc);*/
  
  /* Solve for f_odd in A_odd_trunc.f_odd_trunc = s_odd_trunc */
  gsl_linalg_LU_decomp (A_odd_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_odd_trunc, permute, s_odd_trunc, f_odd_trunc); 
  print_vector(f_odd_trunc);

  /* Set solution vector f_odd = (0, ..., f_odd_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_odd, i+p, gsl_vector_get(f_odd_trunc, i));
  
  /* free memory */
  gsl_matrix_free(A_odd);
  gsl_matrix_free(A_odd_trunc);
  gsl_vector_free(s_odd_trunc);
  gsl_vector_free(f_odd_trunc);
  gsl_permutation_free(permute);

}


/*******************************************************/
/* Solve shell:                                        */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_shell(int L, double alpha, double beta, gsl_vector *s_old, gsl_vector *f)
{
  int N = (*s_old).size;
  int i, j;
  
  gsl_matrix *A = gsl_matrix_alloc (N, N);
  gsl_matrix *x = gsl_matrix_alloc (N, N);
  gsl_matrix *id = gsl_matrix_alloc (N, N);
  gsl_matrix *axplusb = gsl_matrix_calloc (N, N); /* set to zero */
  gsl_matrix *axplusb2 = gsl_matrix_alloc (N, N);

  gsl_vector *s = gsl_vector_alloc (N);
 
  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N);
  
  gsl_vector_set_zero(f);
  
  /* make operator matrix */
  set_A_kernel_even(A, L);
  
  /* multiply s by (alpha x + beta)^2 */
  set_x(x);
  gsl_matrix_scale(x, alpha);
  gsl_matrix_set_identity(id);
  gsl_matrix_scale(id, beta);
  gsl_matrix_add(axplusb, x);
  gsl_matrix_add(axplusb, id);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, axplusb, axplusb,
				 0.0, axplusb2);
  gsl_blas_dgemv(CblasNoTrans,
				 1.0, axplusb2, s_old,
				 0.0, s);

  /* make matrix banded and perform the same operations on source vector */
  /* L_i = \frac{(1+\delta_{0i})L_i - L_{i+2}}{i+1}, for 0 <= i <= N-3   */
  /* L_i = L_i - L_{i+2},                            for 0 <= i <= N-5   */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, 
					 ((1+delta(i, 0))*gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j))/(i+1));
	gsl_vector_set(s, i, 
				   ((1+delta(i, 0))*gsl_vector_get(s, i) - gsl_vector_get(s, i+2))/(i+1));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
  }
  /*printf("A_even banded:\n");
	print_matrix(A);
	print_vector(s);*/
  
  /* Solve for f in A.f = s */
  gsl_linalg_LU_decomp (A, permute, &k); 
  gsl_linalg_LU_solve (A, permute, s, f); 
  print_vector(f);
			   
  /* free memory */
  gsl_matrix_free(A);
  gsl_permutation_free(permute);
}


/*******************************************************/
/* Solve external domain:                              */
/* Take a source vector and moment L and output the    */
/* solution vector                                     */
/*******************************************************/
void solve_ext(int L, double alpha, gsl_vector *s_old, gsl_vector *f)
{
  int N = (*s_old).size;
  int i, j;
  int p; /* order of eigenvalue for homogeneous solution to A */
  
  if(L == 0)
	p = 2;
  else
	p = 3;
  
  gsl_matrix *A = gsl_matrix_alloc (N, N);
  gsl_matrix *xmin1inv = gsl_matrix_alloc (N, N);
  gsl_matrix *xmin1inv2 = gsl_matrix_alloc (N, N);
  gsl_matrix *xmin1inv4 = gsl_matrix_alloc (N, N);
  gsl_matrix *A_trunc = gsl_matrix_alloc (N-p, N-p);
  gsl_vector *s = gsl_vector_alloc (N);
  gsl_vector *s_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_trunc = gsl_vector_alloc (N-p);
  
  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);
  
  gsl_vector_set_zero(f);
  
  /* make operator matrix */
  set_A_ext(A, L);
  
  /* divide s by \alpha^2(x-1)^4 */
  set_xmin1inv(xmin1inv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, xmin1inv, xmin1inv,
				 0.0, xmin1inv2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
				 1.0, xmin1inv2, xmin1inv2,
				 0.0, xmin1inv4);
  gsl_blas_dgemv(CblasNoTrans,
				 1.0, xmin1inv4, s_old,
				 0.0, s);
  gsl_vector_scale(s, 1/(alpha*alpha));
  
  /* make matrix banded and perform the same operations on source vector */
  /* L_i = (1+\delta_{0i})L_i - L_{i+2}, for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, 
					 (1+delta(i, 0))*gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+1, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+1));
  }
  for(i=0; i<=N-5; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
  }
  /*printf("A banded:\n");
  print_matrix(A);
  print_vector(s);*/
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_trunc, i, j, gsl_matrix_get(A, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_trunc, i, gsl_vector_get(s, i));

  /*printf("A banded truncated:\n");
  print_matrix(A_trunc);
  print_vector(s_trunc);*/

  /* Solve for f_even in A_even_trunc.f_even_trunc = s_even_trunc */
  gsl_linalg_LU_decomp (A_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_trunc, permute, s_trunc, f_trunc); 
  /*print_vector(f_trunc);*/

  /* Set solution vector f = (0, ..., f_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f, i+p, gsl_vector_get(f_trunc, i));
  /*print_vector(f);*/
  
  /* free memory */
  gsl_matrix_free(A); 
  gsl_matrix_free(xmin1inv); 
  gsl_matrix_free(xmin1inv2);
  gsl_matrix_free(xmin1inv4);
  gsl_matrix_free(A_trunc);
  gsl_vector_free(s);
  gsl_vector_free(s_trunc);
  gsl_vector_free(f_trunc);
  gsl_permutation_free(permute);
}


/********************************************************************/
/* Find coefficients for homogeneous solution to satisfy continuity */
/* given the particular solution.                                   */
/* Returns homo_coeff = (A0, A1, B1, A2, B2, ..., Az-2, Bz-2, Bz-1) */
/********************************************************************/
void solve_continuity(int L, gsl_vector *boundary, gsl_matrix *particular_m, gsl_vector *homo_coeff)
{
  int i;
  int k;
  int z = (*particular_m).size1; /* number of zones (must be >= 2) */
  int N = (*particular_m).size2;
  
  double fin;
  double dfin;
  double fout;
  double dfout;
  
  gsl_matrix *cont_m = gsl_matrix_calloc (2*z-2, 2*z-2); /* set elements to zero */
  gsl_vector *cont_v = gsl_vector_calloc (2*z-2);

  int luint;
  gsl_permutation *permute = gsl_permutation_alloc (2*z-2);
  
  /* z zones */  
  /* z-1 boundaries between zones */
  /* z-2 is index of last boundary in boundary vector */
  /* 0 is index of boundary between kernal and 1st shell */
  /* 1...z-3 are indices of boundaries between two shells */
  /* z-2 is index of boundary between last shell and external domain */
  
  /* set values for cont_m */
  if(z < 2){ /* error */
	
	printf("There must be at least 2 zones");
	
  }else if(z == 2){ /* there is only a kernel and an external domain */
	
	/* row for continuity */
	gsl_matrix_set(cont_m, 0, 0, pow(gsl_vector_get(boundary, 0), L));
	gsl_matrix_set(cont_m, 0, 1, -pow(gsl_vector_get(boundary, 0), -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 1, 0, L*pow(gsl_vector_get(boundary, 0), L-1));
	gsl_matrix_set(cont_m, 1, 1, (L+1)*pow(gsl_vector_get(boundary, 0), -L-2));
	
  }else{ /* there are shells */
	
	/* internal boundary: */
	/* row for continuity */
	gsl_matrix_set(cont_m, 0, 0, pow(gsl_vector_get(boundary, 0), L));
	gsl_matrix_set(cont_m, 0, 1, -pow(gsl_vector_get(boundary, 0), L));
	gsl_matrix_set(cont_m, 0, 2, -pow(gsl_vector_get(boundary, 0), -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 1, 0, L*pow(gsl_vector_get(boundary, 0), L-1));
	gsl_matrix_set(cont_m, 1, 1, -L*pow(gsl_vector_get(boundary, 0), L-1));
	gsl_matrix_set(cont_m, 1, 2, (L+1)*pow(gsl_vector_get(boundary, 0), -L-2));
	
	/* shell boundaries: */
	
	for(i=1; i<=z-3; i++){
	  /* row for continuity */
	  gsl_matrix_set(cont_m, 2*i, 2*i-1, pow(gsl_vector_get(boundary, i), L));
	  gsl_matrix_set(cont_m, 2*i, 2*i, pow(gsl_vector_get(boundary, i), -L+1));
	  gsl_matrix_set(cont_m, 2*i, 2*i+1, -pow(gsl_vector_get(boundary, i), L));
	  gsl_matrix_set(cont_m, 2*i, 2*i+2, -pow(gsl_vector_get(boundary, i), -L-1));
	  /* row for continuity of 1st derivative */
	  gsl_matrix_set(cont_m, 2*i+1, 2*i-1, L*pow(gsl_vector_get(boundary, i), L-1));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i, -(L+1)*pow(gsl_vector_get(boundary, i), -L-2));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i+1, -L*pow(gsl_vector_get(boundary, i), L-1));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i+2, (L+1)*pow(gsl_vector_get(boundary, i), -L-2));	
	}
	
	/* external boundary: */
	/* row for continuity */
	gsl_matrix_set(cont_m, 2*z-4, 2*z-5, pow(gsl_vector_get(boundary, z-2), L));
	gsl_matrix_set(cont_m, 2*z-4, 2*z-4, pow(gsl_vector_get(boundary, z-2), -L-1));
	gsl_matrix_set(cont_m, 2*z-4, 2*z-3, -pow(gsl_vector_get(boundary, z-2), -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 2*z-3, 2*z-5, L*pow(gsl_vector_get(boundary, z-2), L-1));
	gsl_matrix_set(cont_m, 2*z-3, 2*z-4, -(L+1)*pow(gsl_vector_get(boundary, z-2), -L-2));
	gsl_matrix_set(cont_m, 2*z-3, 2*z-3, (L+1)*pow(gsl_vector_get(boundary, z-2), -L-2));
	
  }
  
  /*print_matrix(cont_m);*/
  
  /* set values for cont_v */
  for(i=0; i<=z-2; i++){
	fin = dfin = fout = dfout = 0;
	for(k=0; k<N; k++){
	  fin += gsl_matrix_get(particular_m, i, k);
	  if(i==0 && !(L%2)){ /* if zone is kernel and L is even */
		dfin += 4*k*k*gsl_matrix_get(particular_m, i, k);
	  } else if(i==0 && L%2){ /* if zone is kernel and L is odd */
		dfin += (2*k+1)*(2*k+1)*gsl_matrix_get(particular_m, i, k);
	  } else { /* not kernel */
		dfin += k*k*gsl_matrix_get(particular_m, i, k);
	  }
	  fout += neg1toi(k)*gsl_matrix_get(particular_m, i+1, k);
	  dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, i+1, k);
	}
	gsl_vector_set(cont_v, 2*i, -fin+fout);
	gsl_vector_set(cont_v, 2*i+1, -dfin+dfout);
  }

  /*print_vector(cont_v);*/
  
  /* Solve for homo_coeff in cont_m.homo_coeff = cont_v */
  gsl_linalg_LU_decomp (cont_m, permute, &luint); 
  gsl_linalg_LU_solve (cont_m, permute, cont_v, homo_coeff); 

  /* free memory */
  gsl_matrix_free(cont_m);
  gsl_vector_free(cont_v);
  gsl_permutation_free(permute);
}


/**********************************************************************/
/* Find the coefficients for the radial part of the Poisson equation. */
/**********************************************************************/
void solve_radial(int L, 
				  int N,
				  double (*source)(double), 
				  double (*source_ext)(double), 
				  gsl_vector *boundary, 
				  gsl_matrix *particular_m, 
				  gsl_vector *homo_coeff)
{
  int nz = (*boundary).size+1; /* number of zones (must be >= 2) */
  int z; /* particular zone being considered */
  int i;
  
  double rlow; /* boundary just less than r */
  double rhigh; /* boundary just greater than r */
  double ulow; /* lower boundary in external domain (u=1/r) */
  double alpha;
  double beta;
  
  gsl_vector *s_kernel = gsl_vector_alloc (N);
  gsl_vector *s_ext = gsl_vector_alloc (N);
  gsl_vector *s_shell = gsl_vector_alloc (N);
  gsl_vector *s_kernel_even = gsl_vector_calloc (2*N-1);
  gsl_vector *s_kernel_odd = gsl_vector_calloc (2*N);
  
  gsl_vector *f_kernel = gsl_vector_alloc (N);
  gsl_vector *f_ext = gsl_vector_alloc (N);
  gsl_vector *f_shell = gsl_vector_alloc (N);
  gsl_vector *f_kernel_even = gsl_vector_calloc (2*N-1);
  gsl_vector *f_kernel_odd = gsl_vector_calloc (2*N);
  
  /* decompose and solve kernel when L is even: */
  if(!(L%2)){ 
    
    rhigh = gsl_vector_get(boundary, 0);
    alpha = rhigh;
    /* go from r in [0, rhigh] to x in [0, 1] domain */
    chebfit(-rhigh, rhigh, s_kernel_even, source);
    /* only use even terms */
    for(i=0; i<N; i++)
      gsl_vector_set(s_kernel, i, gsl_vector_get(s_kernel_even, 2*i)); 
    solve_kernel_even(L, alpha, s_kernel, f_kernel);  
    
    /* decompose and solve kernel when L is odd: */
  } else {

    rhigh = gsl_vector_get(boundary, 0);
    alpha = rhigh;
    /* go from r in [0, rhigh] to x in [0, 1] domain */
    chebfit(-rhigh, rhigh, s_kernel_odd, source);
    /* only use even terms */
    for(i=0; i<N; i++)
      gsl_vector_set(s_kernel, i, gsl_vector_get(s_kernel_odd, 2*i+1)); 
    solve_kernel_even(L, alpha, s_kernel, f_kernel);  
    
  }
  for(i=0; i<N; i++){
    gsl_matrix_set(particular_m, 0, i, gsl_vector_get(f_kernel, i)); 
  }

  /* decompose and solve external domain: */
  ulow = gsl_vector_get(boundary, nz-2);
  alpha = -0.5/ulow;
  /* go from u in [ulow, 0] to x in [-1, 1] domain */
  chebfit(ulow, 0, s_ext, source_ext);
  solve_ext(L, alpha, s_ext, f_ext); 
  for(i=0; i<N; i++){
    gsl_matrix_set(particular_m, nz-1, i, gsl_vector_get(f_ext, i)); 
  }  
  
  /* decompose and solve shell(s) if there are any: */
  if(nz>=3){
    for(z=1; z<=nz-2; z++){
      rlow = gsl_vector_get(boundary, z-1);
      rhigh = gsl_vector_get(boundary, z);
      alpha = 0.5*(rhigh - rlow);
      beta = 0.5*(rhigh + rlow);
      /* go from r in [rlow, rhigh] to x in [-1, 1] domain */
      chebfit(rlow, rhigh, s_shell, source);
      solve_shell(L, alpha, beta, s_shell, f_shell); 
      for(i=0; i<N; i++){
		gsl_matrix_set(particular_m, z, i, gsl_vector_get(f_shell, i)); 
      }
    }
  }  
  
  /* require continuity of function and derivative: */ 
  solve_continuity(L, boundary, particular_m, homo_coeff);  
  
  /* free memory */
  gsl_vector_free(s_kernel);
  gsl_vector_free(s_ext);
  gsl_vector_free(s_shell);
  gsl_vector_free(s_kernel_even);
  gsl_vector_free(s_kernel_odd);
  
  gsl_vector_free(f_kernel);
  gsl_vector_free(f_ext);
  gsl_vector_free(f_shell);
  gsl_vector_free(f_kernel_even);
  gsl_vector_free(f_kernel_odd);
}



/****************************************************************/
/* Evaluate the radial part of the Poisson equation at point r. */
/****************************************************************/
double radialeval(int L, 
				  gsl_vector *boundary, 
				  gsl_matrix *particular_m, 
				  gsl_vector *homo_coeff, 
				  double r)
{
  int nz = (*particular_m).size1; /* number of zones (must be >= 2) */
  int N = (*particular_m).size2; /* number of polynomials */
  int ihigh, ilow, imid;
  int z; /* particular zone r is in */
  int i;
  double a, b; /* coefficients for r^L and r^-(L+1) */
  double rlow; /* boundary just less than r */
  double rhigh; /* boundary just greater than r */

  /* determine what zone r is in (labeled 0, 1, ..., nz-2, nz-1) */
  /* (boundaries between zones are labeled 0, 1, ..., nz-2) */
  if(r <= gsl_vector_get(boundary, 0)) {
	/* r is in kernel */
	z=0;
  } else if(r >= gsl_vector_get(boundary, nz-2)) {
	/* r is in external domain */
	z=nz-1;
  } else {
	/* r is in a shell */
	/* do binary search to find which shell r is in */
	ilow = 0;
	ihigh = z-2; /* z-1 boundaries.  z-2 is index of last boundary. */
	imid = (ihigh+ilow)/2; /* rounds down when half integer */
	while((ihigh-ilow) > 1) {
	  if(r >= gsl_vector_get(boundary, imid)) {
		ilow = imid;
		imid = (ihigh+imid)/2;
	  } else {
      ihigh = imid;
      imid = (ilow+imid)/2;
	  }
	} /* boundaries ilow and ihigh bracket point r */
	z=ihigh;
  }
  
  /* if zone is kernel and L is even: */
  if(z==0 && !(L%2)){ 
	a = gsl_vector_get(homo_coeff, 0);
	rhigh = gsl_vector_get(boundary, 0);
	gsl_vector *coeff = gsl_vector_calloc (2*N-1);
	for(i=0; i<N; i++){
	  gsl_vector_set(coeff, 2*i, gsl_matrix_get(particular_m, 0, i));
	}
	/*print_vector(coeff);*/
	return chebeval(-rhigh, rhigh, coeff, r) + a*pow(r, L);
	/* if zone is kernel and L is odd: */
  } else if(z==0 && L%2){
	a = gsl_vector_get(homo_coeff, 0);
	rhigh = gsl_vector_get(boundary, 0);
	gsl_vector *coeff = gsl_vector_calloc (2*N);
	for(i=0; i<N; i++){
	  gsl_vector_set(coeff, 2*i+1, gsl_matrix_get(particular_m, 0, i));
	}
	return chebeval(-rhigh, rhigh, coeff, r) + a*pow(r, L);
	/* if zone is external domain: */
  } else if(z==nz-1){ 
	b = gsl_vector_get(homo_coeff, 2*nz-3);
	rlow = gsl_vector_get(boundary, nz-2);
	gsl_vector *coeff = gsl_vector_calloc (N);
	for(i=0; i<N; i++){
	  gsl_vector_set(coeff, i, gsl_matrix_get(particular_m, nz-1, i));
	}
	return chebeval(1/rlow, 0, coeff, 1/r) + b*pow(r, -L-1);
	/* zone is shell: */
  } else { 
	a = gsl_vector_get(homo_coeff, 2*z-1);
	b = gsl_vector_get(homo_coeff, 2*z);
	rlow = gsl_vector_get(boundary, nz-2);
	rhigh = gsl_vector_get(boundary, nz-2);
	gsl_vector *coeff = gsl_vector_calloc (N);
	for(i=0; i<N; i++){
	  gsl_vector_set(coeff, i, gsl_matrix_get(particular_m, z-1, i));
	}
	return chebeval(rlow, rhigh, coeff, r) + a*pow(r, L) + b*pow(r, -L-1);
  }
}


/*******************************************************/
/* Going between functions and Chebyshev coefficients. */
/*******************************************************/


/*************************************************************/
/* Takes a function and finds N Chebyshev coefficients       */
/* upto order N-1 using Chebyshev-Gauss-Lobatto quadrature.  */
/*************************************************************/
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

/****************************************************************/
/* Takes N coefficients of Chebyshev polynomials upto order N-1 */
/* and uses Clenshaw's recurrance formula to determine the      */
/* function value at x.                                         */
/****************************************************************/
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
