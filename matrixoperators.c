/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "poisson.h"

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
