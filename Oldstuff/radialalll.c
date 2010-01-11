/* To compile type: gcc -g -lm -lgsl -lgslcblas -Wall -pedantic -ansi radialalll.c */

/********************************************************************/
/* TO DO: 1) chebfit_odd, chebeval_odd functions are not working    */
/*        2) use LAPACK to solve pentadiagonal matrices             */
/********************************************************************/

/* make clist structure
   void clist_set(clist *c, int L, int m, int n, int z, double value);
   double clist_get(clist *c, int L, int m, int n, int z); */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* constants and macros */
#define PI 3.141592653589793
#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
#define cosipiby2(i) ((i)%2 ? 0 : ((i)%4==0 ? 1 : -1)) /* \cos(i\pi/2) */

/* function prototypes */
/*double f_test(double r);*/

double source(double r);
double source_shell(double r);
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
void chebfit_even(double a, gsl_vector *c, double (*func)(double));
void chebfit_odd(double a, gsl_vector *c, double (*func)(double));

double chebeval(double a, double b, gsl_vector *c, double x);
double chebeval_even(double a, gsl_vector *c, double x);
double chebeval_odd(double a, gsl_vector *c, double x);

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
  int N = 8;
  int L = 4;
  
  int nz = 5; /* number of zones (must be >= 2) */
  gsl_vector *homo_coeff = gsl_vector_alloc (2*nz-2);
  gsl_vector *boundary = gsl_vector_alloc (nz-1); /* position of each boundary */
  /*                |Chebyshev coefficients of solution in zone 0   | */
  /*  		        |Chebyshev coefficients of solution in zone 1   | */
  /* particular_m = |                       ...                     | */
  /*                |Chebyshev coefficients of solution in zone nz-1| */
  gsl_matrix *particular_m = gsl_matrix_alloc (nz, N);
   
  double r;
  
  FILE *fp;
  
/*   gsl_vector *test_v = gsl_vector_alloc (8); */
  
/*   chebfit_odd(10.0, test_v, f_test); */
/*   print_vector(test_v); */
/*   for(r = 0.0; r <=10.0; r += 0.1) { */
/* 	printf("%f\t%.18e\n", r, chebeval_odd(10.0, test_v, r)); */
/*   } */
  
/*   gsl_vector_set(boundary, 0, 5.0); */
  
  gsl_vector_set(boundary, 0, 1.0);
  gsl_vector_set(boundary, 1, 2.0);
  gsl_vector_set(boundary, 2, 4.0); 
  gsl_vector_set(boundary, 3, 5.0);

/*   gsl_vector_set(boundary, 0, 1.0); */
/*   gsl_vector_set(boundary, 1, 3.0); */
/*   gsl_vector_set(boundary, 2, 5.0); */
/*   gsl_vector_set(boundary, 3, 10.0); */
/*   gsl_vector_set(boundary, 4, 15.0); */

  solve_radial(L, N, source, source_ext, boundary, particular_m, homo_coeff); 
  
  printf("boundary:\n");
  print_vector(boundary); 
  printf("particular solution matrix:\n");
  print_matrix(particular_m); 
  
  fp=fopen("fofr.txt", "w");
  for(r = 0.0; r <=20.0; r += 0.005) {
    fprintf(fp, "%.18e\t%.18e\n", r, radialeval(L, boundary, particular_m, homo_coeff, r));
	/*printf("%f\t%.18e\n", r, radialeval(L, boundary, particular_m, homo_coeff, r));*/
  }
  
  /* free memory */
  gsl_vector_free(homo_coeff);
  gsl_vector_free(boundary);
  gsl_matrix_free(particular_m);
  
  return 0;
}


/* double f_test(double r) */
/* { */
/*   return -4*r+40*pow(r,3)-96*pow(r,5)+64*pow(r,7); */
/* } */


/*************************************/
/*            Source                 */
/*************************************/

/* double source(double r) */
/* { */
/*   double surface = 5.0; */
/*   return surface - r*r/surface; */
/* } */

/* double source_ext(double u) */
/* { */
/*   double surface = 5.0; */
/*   return pow(surface, 5)*pow(u, 4); */
/* } */

/* double source(double r) */
/* { */
/*   int L = 4; */
/*   double R = 5.0; */
/*   if(r<R){ */
/* 	return 0.5*(L+3)*(L+5)*((L-4)*r*r/pow(R,L+5)-(L-2)/pow(R,L+3)); */
/*   }else{ */
/* 	return 0.0; */
/*   } */
/* } */

/* double source(double r) */
/* { */
/*   double R = 5.0; */
/*   return sin(2.0*PI*r/R); */
/* } */

double source(double r)
{
  int L = 4;
  double R = 5.0;
  return 0.5*(L+3)*(L+5)*((L-4)*r*r/pow(R,L+5)-(L-2)/pow(R,L+3));
}

/* double source(double r) */
/* { */
/*   return 1.0; */
/* } */

double source_shell(double r)
{
  int L = 4;
  double R = 5.0;
  return 0.5*(L+3)*(L+5)*((L-4)*r*r/pow(R,L+5)-(L-2)/pow(R,L+3));
}

/* double source_shell(double r) */
/* { */
/*   return 1.0; */
/* } */

double source_ext(double u)
{
  return 0.0;
}

/* double source(double r) */
/* { */
/*   int L = 51; */
/*   double R = 5.0; */
/*   return 0.5*(L+4)*(L+6)*(-(L-3)*r/pow(R,L+4)+(L-5)*r*r*r/pow(R,L+6)); */
/* } */

/* double source_ext(double u) */
/* { */
/*   return 0; */
/* } */

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
	printf("L is not even in the function set_A_kernel_even.\n");

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

  double fzero; /* f(x=0) */

  gsl_vector_set_zero(f_even);

  /* make operator matrix */
  set_A_kernel_even(A_even, L);
  printf("A_even:\n");
  print_matrix(A_even);
  print_vector(s_even);

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
  printf("A_even banded:\n");
  print_matrix(A_even);
  print_vector(s_even);
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of source vector                      */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_even_trunc, i, j, gsl_matrix_get(A_even, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_even_trunc, i, gsl_vector_get(s_even, i));
  printf("A_even banded truncated:\n");
  print_matrix(A_even_trunc);
  print_vector(s_even_trunc);
  
  /********************************************************/
  /* TO DO: Use LAPACK to find f instead of what's below, */
  /*        to take advantage of the banded matrix        */
  /********************************************************/

  /* Solve for f_even in A_even_trunc.f_even_trunc = s_even_trunc */
  gsl_linalg_LU_decomp (A_even_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_even_trunc, permute, s_even_trunc, f_even_trunc); 
  /*print_vector(f_even_trunc);*/
  printf("f_even_trunc:\n");
  print_vector(f_even_trunc);


  /* Set solution vector f_even = (0, ..., f_even_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_even, i+p, gsl_vector_get(f_even_trunc, i));
  printf("f_even after restoring leading zeros:\n");
  print_vector(f_even);

  /* satisfy boundary conditions */
  /* f(x=r=0)=0 for L>=2 */
  /* df/dr(r=0)=0 is already satisfied for even L */
  if(L>=2){
	fzero = 0.0;
	for(i=0; i<N; i++)
	  fzero += neg1toi(i)*gsl_vector_get(f_even, i);
	gsl_vector_set(f_even, 0, -fzero); /* first element of f was already zero */
  }
  printf("f_even after requiring solution equal 0 at x=0:\n");
  print_vector(f_even);

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

  double dfzero;

  gsl_vector_set_zero(f_odd);

  /* make operator matrix */
  set_A_kernel_odd(A_odd, L);
  printf("A_kernel_odd:\n");
  print_matrix(A_odd);	
  
  /* scale s_kernel_even by alpha^2 */
  gsl_vector_scale(s_odd, alpha*alpha);
  print_vector(s_odd);

  /* make matrix banded and perform the same operations on source vector */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A_odd, i, j,gsl_matrix_get(A_odd, i, j) - gsl_matrix_get(A_odd, i+2, j));
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
  printf("A_kernel_odd banded:\n");
  print_matrix(A_odd);
  print_vector(s_odd);
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_odd_trunc, i, j, gsl_matrix_get(A_odd, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_odd_trunc, i, gsl_vector_get(s_odd, i));
   
  printf("A_kernel_odd banded truncated:\n");
  print_matrix(A_odd_trunc);
  print_vector(s_odd_trunc);
  
  /* Solve for f_odd in A_odd_trunc.f_odd_trunc = s_odd_trunc */
  gsl_linalg_LU_decomp (A_odd_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_odd_trunc, permute, s_odd_trunc, f_odd_trunc); 
  printf("f_odd_trunc:\n");
  print_vector(f_odd_trunc);

  /* Set solution vector f_odd = (0, ..., f_odd_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_odd, i+p, gsl_vector_get(f_odd_trunc, i));
  printf("f_odd after restoring leading zeros:\n");
  print_vector(f_odd);

  /* satisfy boundary conditions */
  /* df/dr(r=0)=0 for L>=3 */
  /* f(x=r=0)=0 is already satisfied for odd L */
  if(L>=3){
	dfzero = 0.0;
	for(i=0; i<N; i++){
	  dfzero += neg1toi(i)*(2*i+1)*gsl_vector_get(f_odd, i);
	  printf("%f\n", dfzero);
	}
	gsl_vector_set(f_odd, 0, -dfzero); /* first element of f was already zero */
  }
  printf("f_odd after making slope 0 at x=0:\n");
  print_vector(f_odd);
  
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
  int p=2; /* order of eigenvalue for homogeneous solution to A */
  
  gsl_matrix *A = gsl_matrix_alloc (N, N);
  gsl_matrix *x = gsl_matrix_alloc (N, N);
  gsl_matrix *id = gsl_matrix_alloc (N, N);
  gsl_matrix *axplusb = gsl_matrix_calloc (N, N); /* set to zero */
  gsl_matrix *axplusb2 = gsl_matrix_alloc (N, N);

  gsl_vector *s = gsl_vector_alloc (N);
 
  gsl_matrix *A_trunc = gsl_matrix_alloc (N-p, N-p);
 
  gsl_vector *hello = gsl_vector_calloc (N-p);

  gsl_vector *s_trunc = gsl_vector_alloc (N-p);
  gsl_vector *f_trunc = gsl_vector_alloc (N-p);
  
  int k;
  gsl_permutation *permute = gsl_permutation_alloc (N-p);
  
  gsl_vector_set_zero(f);

  printf("alpha=%f\tbeta=%f\n", alpha, beta);
  
  /* make operator matrix */
  set_A_shell(A, L, alpha, beta);
  /*set_A_shell(A, 7, 3.345, 73.2345);*/
  /*gsl_matrix_set_identity(A);*/
  printf("A shell:\n");
  print_matrix(A);
  
  printf("original source before multiplying by (alpha x + beta)^2:\n"); 
  print_vector(s_old);
  
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
  printf("source after multiplying by (alpha x + beta)^2:\n"); 
  print_vector(s);

  /* make matrix banded and perform the same operations on source vector */
  /* L_i = \frac{(1+\delta_{0i})L_i - L_{i+2}}{i+1}, for 0 <= i <= N-3   */
  /* L_i = L_i - L_{i+2},                            for 0 <= i <= N-5   */
  /* for(i=0; i<=N-3; i++) { */
/* 	for(j=0; j<N; j++) */
/* 	  gsl_matrix_set(A, i, j, */
/* 					 ((1+delta(i, 0))*gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j))/(i+1)); */
/* 	gsl_vector_set(s, i, */
/* 				   ((1+delta(i, 0))*gsl_vector_get(s, i) - gsl_vector_get(s, i+2))/(i+1)); */
/*   } */
/*   for(i=0; i<=N-5; i++) { */
/* 	for(j=0; j<N; j++) */
/* 	  gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j)); */
/* 	gsl_vector_set(s, i, gsl_vector_get(s, i) - gsl_vector_get(s, i+2)); */
/*   } */
/*   printf("A shell banded:\n"); */
/*   print_matrix(A); */
/*   print_vector(s); */
  

  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_trunc, i, j, gsl_matrix_get(A, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_trunc, i, gsl_vector_get(s, i));
  

  /*gsl_matrix_set_identity(A_trunc);*/


  printf("A shell banded truncated:\n");
  print_matrix(A_trunc);
  print_vector(s_trunc);
  
  
  
  /* Solve for f_trunc in A_trunc.f_trunc = s_trunc */
  gsl_linalg_LU_decomp (A_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_trunc, permute, s_trunc, f_trunc); 
  printf("f_trunc for shell:\n");
  print_vector(f_trunc);
  

  
  gsl_blas_dgemv(CblasNoTrans,
				 1.0, A_trunc, f_trunc,
				 0.0, hello);
  printf("A_trunc.f_trunc should equal s_trunc:\n");
  print_vector(hello);
  



  /* Set solution vector f = (0, ..., f_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f, i+p, gsl_vector_get(f_trunc, i));
  printf("f after restoring leading zeros:\n");
  print_vector(f);
  

  /* free memory */
  gsl_matrix_free(A);
  gsl_matrix_free(x);
  gsl_matrix_free(id);
  gsl_matrix_free(axplusb);
  gsl_matrix_free(axplusb2);
  gsl_matrix_free(A_trunc);
  
  gsl_vector_free(s);
  gsl_vector_free(s_trunc);
  gsl_vector_free(f_trunc);
    
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
  
  double finf;
  double dfinf;

  gsl_vector_set_zero(f);
  
  /* make operator matrix */
  set_A_ext(A, L);
  printf("A_ext:\n");
  print_matrix(A);

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
  
  printf("LHS of equation for external domain:\n");
  print_vector(s);

  /* make matrix banded and perform the same operations on source vector */
  /* L_i = (1+\delta_{0i})L_i - L_{i+2}, for 0 <= i <= N-3               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+1},                for 0 <= i <= N-5               */
  /* L_i = L_i - L_{i+2},                for 0 <= i <= N-5               */
  for(i=0; i<=N-3; i++) {
	for(j=0; j<N; j++)
	  gsl_matrix_set(A, i, j, 
					 (1+delta(i, 0))*gsl_matrix_get(A, i, j) - gsl_matrix_get(A, i+2, j));
	gsl_vector_set(s, i, (1+delta(i, 0))*gsl_vector_get(s, i) - gsl_vector_get(s, i+2));
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
  printf("A banded:\n");
  print_matrix(A);
  print_vector(s);
  
  /* truncate operator matrix and source vector            */
  /* by removing first p columns and last p rows of matrix */
  /* and last p rows of vector                             */
  for(i=0; i<N-p; i++) {
	for(j=0; j<N-p; j++)
	  gsl_matrix_set(A_trunc, i, j, gsl_matrix_get(A, i, j+p));
  }
  for(i=0; i<N-p; i++)
	gsl_vector_set(s_trunc, i, gsl_vector_get(s, i));

  printf("A banded truncated:\n");
  print_matrix(A_trunc);
  print_vector(s_trunc);
  
  /* Solve for f_even in A_even_trunc.f_even_trunc = s_even_trunc */
  gsl_linalg_LU_decomp (A_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_trunc, permute, s_trunc, f_trunc); 
  print_vector(f_trunc);
  
  /* Set solution vector f = (0, ..., f_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f, i+p, gsl_vector_get(f_trunc, i));
  print_vector(f);
 
  /* satisfy boundary conditions */
  /* f_particular(x=1) + \alpha*T_0(x=1) = 0 */
  finf = 0;
  for(i=2; i<N; i++)
	finf += gsl_vector_get(f, i);
  gsl_vector_set(f, 0, -finf); /* first element of f was already zero */
  /* set df/dx=0 for L>=1 by changing T_0 term */
  if(L>=1){
	dfinf = 0;
	for(i=2; i<N; i++)
	  dfinf += i*i*gsl_vector_get(f, i);
	gsl_vector_set(f, 1, -dfinf);
  }
    
  printf("external solution before homogeneous part is added\n");
  printf("but after T_0, T_1 terms are accounted for:\n");
  print_vector(f);
  
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
void solve_continuity(int L, 
					  gsl_vector *boundary, 
					  gsl_matrix *particular_m, 
					  gsl_vector *homo_coeff)
{
  int i;
  int k;
  int z = (*particular_m).size1; /* number of zones (must be >= 2) */
  int N = (*particular_m).size2;
  
  double R;
  double alpha;
  double fin;
  double dfin;
  double fout;
  double dfout;
  
  gsl_matrix *cont_m = gsl_matrix_calloc (2*z-2, 2*z-2); /* set elements to zero */
  gsl_vector *cont_v = gsl_vector_calloc (2*z-2); /* set elements to zero */

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
	R = gsl_vector_get(boundary, 0);
	/* row for continuity */
	gsl_matrix_set(cont_m, 0, 0, pow(R, L));
	gsl_matrix_set(cont_m, 0, 1, -pow(R, -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 1, 0, L*pow(R, L-1));
	gsl_matrix_set(cont_m, 1, 1, (L+1)*pow(R, -L-2));
	
  }else{ /* there are shells */
	
	/* internal boundary: */
	R = gsl_vector_get(boundary, 0);
	/* row for continuity */
	gsl_matrix_set(cont_m, 0, 0, pow(R, L));
	gsl_matrix_set(cont_m, 0, 1, -pow(R, L));
	gsl_matrix_set(cont_m, 0, 2, -pow(R, -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 1, 0, L*pow(R, L-1));
	gsl_matrix_set(cont_m, 1, 1, -L*pow(R, L-1));
	gsl_matrix_set(cont_m, 1, 2, (L+1)*pow(R, -L-2));
	
	/* shell boundaries: */
	
	for(i=1; i<=z-3; i++){
	  R = gsl_vector_get(boundary, i);
	  /* row for continuity */
	  gsl_matrix_set(cont_m, 2*i, 2*i-1, pow(R, L));
	  gsl_matrix_set(cont_m, 2*i, 2*i, pow(R, -L-1));
	  gsl_matrix_set(cont_m, 2*i, 2*i+1, -pow(R, L));
	  gsl_matrix_set(cont_m, 2*i, 2*i+2, -pow(R, -L-1));
	  /* row for continuity of 1st derivative */
	  gsl_matrix_set(cont_m, 2*i+1, 2*i-1, L*pow(R, L-1));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i, -(L+1)*pow(R, -L-2));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i+1, -L*pow(R, L-1));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i+2, (L+1)*pow(R, -L-2));	
	}
	
	/* external boundary: */
	R = gsl_vector_get(boundary, z-2);
	/* row for continuity */
	gsl_matrix_set(cont_m, 2*z-4, 2*z-5, pow(R, L));
	gsl_matrix_set(cont_m, 2*z-4, 2*z-4, pow(R, -L-1));
	gsl_matrix_set(cont_m, 2*z-4, 2*z-3, -pow(R, -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 2*z-3, 2*z-5, L*pow(R, L-1));
	gsl_matrix_set(cont_m, 2*z-3, 2*z-4, -(L+1)*pow(R, -L-2));
	gsl_matrix_set(cont_m, 2*z-3, 2*z-3, (L+1)*pow(R, -L-2));
	
  }
  
  printf("continuity matrix:\n");
  print_matrix(cont_m);
  
  /* set values for cont_v */
  if(z == 2){ /* there is only a kernel and an external domain */
	
	fin = dfin = fout = dfout = 0; 
	R = gsl_vector_get(boundary, 0);
	/* function and derivative on inside: */
	alpha = R;
	if(!(L%2)){ /* L is even */
	  for(k=0; k<N; k++){	
		fin += gsl_matrix_get(particular_m, 0, k);
		dfin += 4*k*k*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	} else { /* L is odd */
	  for(k=0; k<N; k++){
		fin += gsl_matrix_get(particular_m, 0, k);		
		dfin += (2*k+1)*(2*k+1)*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	}
	/* function and derivative on outside: */	
	alpha = -0.5/R;
	for(k=0; k<N; k++){	
	  fout += neg1toi(k)*gsl_matrix_get(particular_m, 1, k);
	  dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, 1, k);
	}
	dfout /= -(alpha*R*R);
	
	gsl_vector_set(cont_v, 0, -fin+fout);
	gsl_vector_set(cont_v, 1, -dfin+dfout);
	
  }else{ /* there are shells */

	/* interface between kernel and first shell */
	fin = dfin = fout = dfout = 0; 
	/* function and derivative on inside: */
	R = gsl_vector_get(boundary, 0);
	alpha = R;
	if(!(L%2)){ /* L is even */
	  for(k=0; k<N; k++){	
		fin += gsl_matrix_get(particular_m, 0, k);
		dfin += 4*k*k*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	} else { /* L is odd */
	  for(k=0; k<N; k++){
		fin += gsl_matrix_get(particular_m, 0, k);		
		dfin += (2*k+1)*(2*k+1)*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	}
	/* function and derivative on outside: */	
	alpha = 0.5*(gsl_vector_get(boundary, 1)-R);
	for(k=0; k<N; k++){	
	  fout += neg1toi(k)*gsl_matrix_get(particular_m, 1, k);
	  dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, 1, k);
	}
	dfout /= alpha;
	
	gsl_vector_set(cont_v, 0, -fin+fout);
	gsl_vector_set(cont_v, 1, -dfin+dfout);
	
	/* interface between shells */
	/* (skipped if z<=3) */
	for(i=1; i<=z-3; i++){
	  fin = dfin = fout = dfout = 0; 
	  R = gsl_vector_get(boundary, i);
	  /* function and derivative on inside: */
	  alpha = 0.5*(R-gsl_vector_get(boundary, i-1));
	  for(k=0; k<N; k++){	
		fin += gsl_matrix_get(particular_m, i, k);
		dfin += k*k*gsl_matrix_get(particular_m, i, k);
	  }
	  dfin /= alpha;
	  /* function and derivative on outside: */	
	  alpha = 0.5*(gsl_vector_get(boundary, i+1)-R);
	  for(k=0; k<N; k++){	
		fout += neg1toi(k)*gsl_matrix_get(particular_m, i+1, k);
		dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, i+1, k);
	  }
	  dfout /= alpha;
	  
	  gsl_vector_set(cont_v, 2*i, -fin+fout);
	  gsl_vector_set(cont_v, 2*i+1, -dfin+dfout);
	}
	
	/* interface between last shell and external domain */
	fin = dfin = fout = dfout = 0; 
	R = gsl_vector_get(boundary, z-2);
	/* function and derivative on inside: */
	alpha = 0.5*(R-gsl_vector_get(boundary, z-3));
	for(k=0; k<N; k++){	
	  fin += gsl_matrix_get(particular_m, z-2, k);
	  dfin += k*k*gsl_matrix_get(particular_m, z-2, k);
	}
	dfin /= alpha;
	/* function and derivative on outside: */	
	alpha = -0.5/R;
	for(k=0; k<N; k++){	
	  fout += neg1toi(k)*gsl_matrix_get(particular_m, z-1, k);
	  dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, z-1, k);
	}
	dfout /= -1/(alpha*R*R);
	
	gsl_vector_set(cont_v, 2*z-4, -fin+fout);
	gsl_vector_set(cont_v, 2*z-3, -dfin+dfout);
  }
  
  printf("vector for continuity equation:\n");
  print_vector(cont_v);
  
  /* Solve for homo_coeff in cont_m.homo_coeff = cont_v */
  gsl_linalg_LU_decomp (cont_m, permute, &luint); 
  gsl_linalg_LU_solve (cont_m, permute, cont_v, homo_coeff); 
  
  printf("solution to continuity equation:\n");
  print_vector(homo_coeff);

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
  double alpha;
  double beta;
  
  gsl_vector *s_kernel = gsl_vector_alloc (N);
  gsl_vector *s_ext = gsl_vector_alloc (N);
  gsl_vector *s_shell = gsl_vector_alloc (N);
  gsl_vector *s_kernel_odd = gsl_vector_calloc (2*N);
  
  gsl_vector *f_kernel = gsl_vector_alloc (N);
  gsl_vector *f_ext = gsl_vector_alloc (N);
  gsl_vector *f_shell = gsl_vector_alloc (N);
  gsl_vector *f_kernel_odd = gsl_vector_calloc (2*N);
  
  /* decompose and solve kernel when L is even: */
  if(!(L%2)){ 
    
    rhigh = gsl_vector_get(boundary, 0);
    alpha = rhigh;
    /* go from r in [0, rhigh] to x in [0, 1] domain */
    chebfit_even(rhigh, s_kernel, source);
    printf("s_kernel:\n");
	print_vector(s_kernel);
    solve_kernel_even(L, alpha, s_kernel, f_kernel);  
	printf("f_kernel:\n");
	print_vector(f_kernel);
    /* decompose and solve kernel when L is odd: */
  } else {

    rhigh = gsl_vector_get(boundary, 0);
    alpha = rhigh;
    /* go from r in [0, rhigh] to x in [0, 1] domain */
    chebfit(-rhigh, rhigh, s_kernel_odd, source);
	printf("original s_kernel_odd:\n");
	print_vector(s_kernel_odd);
    /* only use even terms */
    for(i=0; i<N; i++)
      gsl_vector_set(s_kernel, i, gsl_vector_get(s_kernel_odd, 2*i+1)); 
	printf("new s_kernel:\n");
	print_vector(s_kernel);
    solve_kernel_odd(L, alpha, s_kernel, f_kernel);  
  }
  
  for(i=0; i<N; i++){
    gsl_matrix_set(particular_m, 0, i, gsl_vector_get(f_kernel, i)); 
  }

  /* decompose and solve external domain: */
  rlow = gsl_vector_get(boundary, nz-2);
  alpha = -0.5/rlow;
  /* go from u in [ulow, 0] to x in [-1, 1] domain */
  chebfit(1/rlow, 0, s_ext, source_ext);
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
	  printf("for shell z=%d, rlow=%g, rhigh=%g, alpha=%g, beta=%g\n", z, rlow, rhigh, alpha, beta);
      /* go from r in [rlow, rhigh] to x in [-1, 1] domain */
      chebfit(rlow, rhigh, s_shell, source_shell);
	  print_vector(s_shell);
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
  gsl_vector_free(s_kernel_odd);
  
  gsl_vector_free(f_kernel);
  gsl_vector_free(f_ext);
  gsl_vector_free(f_shell);
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
  int zhigh, zlow, zmid;
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
	zlow = 0; /* index of zone for first shell */
	zhigh = nz-2; /* index of zone for last shell */
	zmid = (zhigh+zlow)/2; /* rounds down when half integer */
	while((zhigh-zlow) > 1) {
	  if(r >= gsl_vector_get(boundary, zmid)) {
		zlow = zmid;
		zmid = (zhigh+zmid)/2;
	  } else {
      zhigh = zmid;
      zmid = (zlow+zmid)/2;
	  }
	} /* boundaries zlow and zhigh bracket point r */
	z=zhigh;
  }
  
  /* if zone is kernel and L is even: */
  if(z==0 && !(L%2)){ 
	a = gsl_vector_get(homo_coeff, 0);
	rhigh = gsl_vector_get(boundary, 0);
	gsl_vector *coeff = gsl_vector_calloc (N);
	for(i=0; i<N; i++){
	  gsl_vector_set(coeff, i, gsl_matrix_get(particular_m, 0, i));
	}
	return chebeval_even(rhigh, coeff, r) + a*pow(r, L);
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
	rlow = gsl_vector_get(boundary, z-1);
	rhigh = gsl_vector_get(boundary, z);
	gsl_vector *coeff = gsl_vector_calloc (N);
	for(i=0; i<N; i++){
	  gsl_vector_set(coeff, i, gsl_matrix_get(particular_m, z, i));
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
    sum=0.0;
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


void chebfit_even(double a, gsl_vector *c, double (*func)(double))
{
  int N;             /* number of terms in Chebyshev expansion */
  int k;             /* kth coefficient of Chebyshev expansion */
  int j;             /* jth term in summation for coefficient c_k */
  double y_j;        /* Gauss-Lobatto quadrature point */
  double sum;        /* sum of terms in summation for c_k */
  double fac;        /* overall factor in front of summation for c_k */
  
  N = (*c).size;
  
  /*Find each of N coefficients c_k.*/
  for(k=0; k<N; k++){
    sum=0;
    /*Evaluate jth term in summation for kth coefficient.*/
    for(j=0; j<N; j++){
      /* Rescale range to be y in [-1, 1] */
	  /* where r=a*sqrt(0.5(y+1)) is in [0, a]. */
      y_j = cos(PI*j/(N-1));
      sum += (*func)(a*sqrt(0.5*(y_j+1)))*cos(PI*k*j/(N-1))/
		(1.0+delta(0, j)+delta(N-1, j));
    }
    fac = 2.0/((N-1)*(1.0+delta(0, k)+delta(N-1, k)));
    gsl_vector_set(c, k, fac*sum);
  }
}


void chebfit_odd(double a, gsl_vector *c, double (*func)(double))
{
  int N;             /* number of terms in Chebyshev expansion */
  int k;             /* kth coefficient of Chebyshev expansion */
  int j;             /* jth term in summation for coefficient c_k */
  double y_j;        /* Gauss-Lobatto quadrature point */
  double sum;        /* sum of terms in summation for c_k */
  double fac;        /* overall factor in front of summation for c_k */
  
  N = (*c).size;
  
  /*Find each of N coefficients c_k.*/
  for(k=0; k<N; k++){
    sum=0;
    /*Evaluate jth term in summation for kth coefficient.*/
    for(j=0; j<N; j++){
      /* Rescale range to be y in [-1, 1] */
	  /* where r=a*sqrt(0.5(y+1)) is in [0, a]. */
      y_j = cos(PI*j/(N-1));
	  if(j<N-1){
		sum += (*func)(a*sqrt(0.5*(y_j+1)))*cos(PI*k*j/(N-1))/
		  ((1.0+delta(0, j)+delta(N-1, j))*a*sqrt(0.5*(y_j+1)));	
	  }else{ 
		sum += (*func)(a*sqrt(0.5*(y_j+1)))*cos(PI*k*j/(N-1))/
		  ((1.0+delta(0, j)+delta(N-1, j))*a*sqrt(0.5*(y_j+1)));
	  }
	  printf("%f\t%f\n", y_j, sqrt(0.5*(y_j+1)));
    }
	
    fac = 2.0/((N-1)*(1.0+delta(0, k)+delta(N-1, k)));
    gsl_vector_set(c, k, fac*sum);
	/* c now represents the coefficients (which are even) for func(x)/x. */	
  }
  /* use recursion relation to multiply series by x: */
  for(k=0; k<=N-2; k++){
	gsl_vector_set(c,k,0.5*((1+delta(0,k))*gsl_vector_get(c,k)+gsl_vector_get(c,k+1)));
  }
  gsl_vector_set(c,N-1,0.5*gsl_vector_get(c,N-1));
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


double chebeval_even(double a, gsl_vector *c, double x)
{ 
  double t; /* change of variable */ 
  double u; /* another change of variable */
  int N; /* number of terms (T_0, T_2, T_4, T_{2N-2}) */
  int k;
  double y_k;
  double y_kplus1; 
  double y_kplus2;  

  N = (*c).size;

  if (x*(x-a) > 0.0)
    printf("x is not between 0 and a in chebeval_even");
  
  /*shift range from [-a, a] to [-1, 1]*/
  t = x/a;
  /*shift range again from [0, 1] to [-1, 1]*/
  u = 2.0*t*t-1;
  /*Begin Clenshaw's recurrance formula.*/
  y_k = 0.0;
  y_kplus1 = 0.0; 
  y_kplus2 = 0.0;
  for(k=N-1; k>=1; k--){
    y_kplus2 = y_kplus1;
    y_kplus1 = y_k;
    y_k = 2*u*y_kplus1 - y_kplus2 + gsl_vector_get(c, k);
  }

  /*k is now 1*/
  return y_k*u - y_kplus1 + gsl_vector_get(c, 0);
}


double chebeval_odd(double a, gsl_vector *c, double x)
{ 
  double t; /* change of variable */ 
  double u; /* another change of variable */
  int N; /* number of terms (T_1, T_3, T_5, T_{2N-1}) */
  int k;
  double y_k;
  double y_kplus1; 
  double y_kplus2;  

  N = (*c).size;

  if (x*(x-a) > 0.0)
    printf("x is not between 0 and a in chebeval_even");
  
  /*shift range from [-a, a] to [-1, 1]*/
  t = x/a;
  /*shift range again from [0, 1] to [-1, 1]*/
  u = 2.0*t*t-1;
  /*Begin Clenshaw's recurrance formula.*/
  y_k = 0.0;
  y_kplus1 = 0.0; 
  y_kplus2 = 0.0;
  for(k=N-1; k>=1; k--){
    y_kplus2 = y_kplus1;
    y_kplus1 = y_k;
    y_k = 2*u*y_kplus1 - y_kplus2 + gsl_vector_get(c, k);
  }

  /*k is now 1*/
  return x*(y_k*(2*u-1) - y_kplus1 + gsl_vector_get(c, 0));
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
    printf ("%8g\t", gsl_vector_get (v, i));
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
      printf ("%8g\t", gsl_matrix_get (m, i, j));
    printf ("\n");
  } 
  printf("\n");
}
