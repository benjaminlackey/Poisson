/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "poisson.h"


/* void solve_poisson(coeff *f_znlm, coeff *s_znlm, bound_coeff *b_zlm); */

/* void solve_kernel_even(int L, double alpha, gsl_vector *s_even, gsl_vector *f_even); */
/* void solve_kernel_odd(int L, double alpha, gsl_vector *s_odd, gsl_vector *f_odd); */
/* void solve_shell(int L, double alpha, double beta, gsl_vector *s_old, gsl_vector *f); */
/* void solve_ext(int L, double alpha, gsl_vector *s_old, gsl_vector *f); */

/* void solve_continuity(int L,  */
/* 					  gsl_vector *boundary,  */
/* 					  gsl_matrix *particular_m,  */
/* 					  gsl_vector *homo_coeff); */

/* void chebfit(double a, double b, gsl_vector *c, double (*func)(double)); */
/* void chebfit_even(double a, gsl_vector *c, double (*func)(double)); */
/* void chebfit_odd(double a, gsl_vector *c, double (*func)(double)); */

/* double chebeval(double a, double b, gsl_vector *c, double x); */
/* double chebeval_even(double a, gsl_vector *c, double x); */
/* double chebeval_odd(double a, gsl_vector *c, double x); */

/* double homogeneous(int L, double a, double b, double x); */

/* void solve_radial(int L,  */
/* 				  gsl_matrix *fieldlm_zn,  */
/* 				  gsl_matrix *sourcelm_zn,  */
/* 				  gsl_vector *boundlm_z);  */

/* void chebfit_homo_kernel_even(gsl_vector *homo_v, int L, double a, double alpha); */
/* void chebfit_homo_kernel_odd(gsl_vector *homo_v, int L, double a, double alpha); */
/* void chebfit_homo_shell(gsl_vector *homo_v, int L, double a, double b, double alpha, double beta); */
/* void chebfit_homo_ext(gsl_vector *homo_v, int L, double b, double alpha); */



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
  
  gsl_matrix *A_even;
  gsl_matrix *A_even_trunc;

  gsl_vector *s_even_trunc;
  gsl_vector *f_even_trunc;

  int k;
  gsl_permutation *permute;

  double fzero; /* f(x=0) */

  
  if(L%2) {
	printf("L is not even in solve_kernel_even.");
	p = 2;
  } else if(L == 0)
	p = 1;
  else
	p = 2;
  
  
  gsl_vector_set_zero(f_even);
  
  A_even = gsl_matrix_alloc (N, N);
  A_even_trunc = gsl_matrix_alloc (N-p, N-p);
  
  s_even_trunc = gsl_vector_alloc (N-p);
  f_even_trunc = gsl_vector_alloc (N-p);
  
  permute = gsl_permutation_alloc (N-p);
  
  
  /* make operator matrix */
  set_A_kernel_even(A_even, L);
  /*printf("A_even:\n");
  print_matrix(A_even);
  print_vector(s_even);*/

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
  /*printf("f_even_trunc:\n");
	print_vector(f_even_trunc);*/


  /* Set solution vector f_even = (0, ..., f_even_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_even, i+p, gsl_vector_get(f_even_trunc, i));
  /*printf("f_even after restoring leading zeros:\n");
	print_vector(f_even);*/

  /* satisfy boundary conditions */
  /* f(x=r=0)=0 for L>=2 */
  /* df/dr(r=0)=0 is already satisfied for even L */
  if(L>=2){
	fzero = 0.0;
	for(i=0; i<N; i++)
	  fzero += neg1toi(i)*gsl_vector_get(f_even, i);
	gsl_vector_set(f_even, 0, -fzero); /* first element of f was already zero */
  }
  /*printf("f_even after requiring solution equal 0 at x=0:\n");
	print_vector(f_even);*/

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
  
  gsl_matrix *A_odd;
  gsl_matrix *A_odd_trunc;
  
  gsl_vector *s_odd_trunc;
  gsl_vector *f_odd_trunc;
  
  int k;
  gsl_permutation *permute;
  
  double dfzero;
  
  if(!(L%2)) {
	printf("L is not odd in solve_kernel_odd.");
	p = 2;
  } else if(L == 1)
	p = 1;
  else
	p = 2;
  
  gsl_vector_set_zero(f_odd);
  A_odd = gsl_matrix_alloc (N, N);
  A_odd_trunc = gsl_matrix_alloc (N-p, N-p);
  
  s_odd_trunc = gsl_vector_alloc (N-p);
  f_odd_trunc = gsl_vector_alloc (N-p);
  
  permute = gsl_permutation_alloc (N-p);
  
  /* make operator matrix */
  set_A_kernel_odd(A_odd, L);
  /*printf("A_kernel_odd:\n");
	print_matrix(A_odd);	*/
  
  /* scale s_kernel_even by alpha^2 */
  gsl_vector_scale(s_odd, alpha*alpha);
  /*print_vector(s_odd);*/

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
   
  /*printf("A_kernel_odd banded truncated:\n");
  print_matrix(A_odd_trunc);
  print_vector(s_odd_trunc);*/
  
  /* Solve for f_odd in A_odd_trunc.f_odd_trunc = s_odd_trunc */
  gsl_linalg_LU_decomp (A_odd_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_odd_trunc, permute, s_odd_trunc, f_odd_trunc); 
  /*printf("f_odd_trunc:\n");
	print_vector(f_odd_trunc);*/

  /* Set solution vector f_odd = (0, ..., f_odd_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f_odd, i+p, gsl_vector_get(f_odd_trunc, i));
  /*printf("f_odd after restoring leading zeros:\n");
	print_vector(f_odd);*/

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
  /*printf("f_odd after making slope 0 at x=0:\n");
	print_vector(f_odd);*/
  
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

  /*printf("alpha=%f\tbeta=%f\n", alpha, beta);*/
  
  /* make operator matrix */
  set_A_shell(A, L, alpha, beta);
  /*set_A_shell(A, 7, 3.345, 73.2345);*/
  /*gsl_matrix_set_identity(A);*/
  /*printf("A shell:\n");
	print_matrix(A);*/
  
  /*printf("original source before multiplying by (alpha x + beta)^2:\n"); 
	print_vector(s_old);*/
  
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
  /*printf("source after multiplying by (alpha x + beta)^2:\n"); 
	print_vector(s);*/

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


  /*printf("A shell banded truncated:\n");
  print_matrix(A_trunc);
  print_vector(s_trunc);*/
  
  
  
  /* Solve for f_trunc in A_trunc.f_trunc = s_trunc */
  gsl_linalg_LU_decomp (A_trunc, permute, &k); 
  gsl_linalg_LU_solve (A_trunc, permute, s_trunc, f_trunc); 
  /*printf("f_trunc for shell:\n");
	print_vector(f_trunc);*/
  

  
  gsl_blas_dgemv(CblasNoTrans,
				 1.0, A_trunc, f_trunc,
				 0.0, hello);
  /*printf("A_trunc.f_trunc should equal s_trunc:\n");
	print_vector(hello);*/
  



  /* Set solution vector f = (0, ..., f_trunc) (p zeros) */
  for(i=0; i<N-p; i++)
	gsl_vector_set(f, i+p, gsl_vector_get(f_trunc, i));
  /*printf("f after restoring leading zeros:\n");
	print_vector(f);*/
  

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
  
  gsl_matrix *A;
  gsl_matrix *xmin1inv;
  gsl_matrix *xmin1inv2;
  gsl_matrix *xmin1inv4;
  gsl_matrix *A_trunc;
  gsl_vector *s;
  gsl_vector *s_trunc;
  gsl_vector *f_trunc;
  
  int k;
  gsl_permutation *permute;
  
  double finf;
  double dfinf;

  if(L == 0)
	p = 2;
  else
	p = 3;
  
  gsl_vector_set_zero(f);
  
  A = gsl_matrix_alloc (N, N);
  xmin1inv = gsl_matrix_alloc (N, N);
  xmin1inv2 = gsl_matrix_alloc (N, N);
  xmin1inv4 = gsl_matrix_alloc (N, N);
  A_trunc = gsl_matrix_alloc (N-p, N-p);
  s = gsl_vector_alloc (N);
  s_trunc = gsl_vector_alloc (N-p);
  f_trunc = gsl_vector_alloc (N-p);
  permute = gsl_permutation_alloc (N-p);
  
  /* make operator matrix */
  set_A_ext(A, L);
  /*printf("A_ext:\n");
	print_matrix(A);*/
  
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
  
  /*printf("LHS of equation for external domain:\n");
	print_vector(s);*/

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
    
  /*printf("external solution before homogeneous part is added\n");
  printf("but after T_0, T_1 terms are accounted for:\n");
  print_vector(f);*/
  
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
  
  /*printf("continuity matrix:\n");
	print_matrix(cont_m);*/
  
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
  
  /*printf("vector for continuity equation:\n");
	print_vector(cont_v);*/
  
  /* Solve for homo_coeff in cont_m.homo_coeff = cont_v */
  gsl_linalg_LU_decomp (cont_m, permute, &luint); 
  gsl_linalg_LU_solve (cont_m, permute, cont_v, homo_coeff); 
  
  /*printf("solution to continuity equation:\n");
	print_vector(homo_coeff);*/

  /* free memory */
  gsl_matrix_free(cont_m);
  gsl_vector_free(cont_v);
  gsl_permutation_free(permute);
}


/**********************************************************************/
/* Find the coefficients for the radial part of the Poisson equation. */
/**********************************************************************/
void solve_radial(int L, 
				  gsl_matrix *fieldlm_zn, 
				  gsl_matrix *sourcelm_zn, 
				  gsl_vector *boundlm_z)
{
  int nz = sourcelm_zn->size1; /* number of zones (must be >= 2) */
  int nn = sourcelm_zn->size2; /* number of Chebyshev polynomials */
  int z; /* particular zone being considered */
  int i;
  double rlow; /* boundary just less than r */
  double rhigh; /* boundary just greater than r */
  double alpha;
  double beta;
  double a;
  double b;
  gsl_vector *s_kernel = gsl_vector_alloc (nn);
  gsl_vector *s_ext = gsl_vector_alloc (nn);
  gsl_vector *s_shell = gsl_vector_alloc (nn);
  gsl_vector *f_kernel = gsl_vector_alloc (nn);
  gsl_vector *f_ext = gsl_vector_alloc (nn);
  gsl_vector *f_shell = gsl_vector_alloc (nn);
  gsl_vector *homo_coeff = gsl_vector_alloc(2*nz-2);
  gsl_vector *homo_v = gsl_vector_alloc(nn);
  gsl_matrix *particularlm_zn = gsl_matrix_alloc(nz, nn);
  
  /* decompose and solve kernel */
  alpha = gsl_vector_get(boundlm_z, 0);
  for(i=0; i<nn; i++)
	gsl_vector_set(s_kernel, i, gsl_matrix_get(sourcelm_zn, 0, i));
  if(!(L%2)) /* L is even: */
	solve_kernel_even(L, alpha, s_kernel, f_kernel); 
  else /* L is odd: */
	solve_kernel_odd(L, alpha, s_kernel, f_kernel);  
  
  for(i=0; i<nn; i++){
    gsl_matrix_set(particularlm_zn, 0, i, gsl_vector_get(f_kernel, i)); 
  }
  
  /* decompose and solve external domain: */
  rlow = gsl_vector_get(boundlm_z, nz-2);
  alpha = -0.5/rlow;
  for(i=0; i<nn; i++)
	gsl_vector_set(s_ext, i, gsl_matrix_get(sourcelm_zn, nz-1, i));
  solve_ext(L, alpha, s_ext, f_ext); 
  
  for(i=0; i<nn; i++){
    gsl_matrix_set(particularlm_zn, nz-1, i, gsl_vector_get(f_ext, i)); 
  }  
  
  /* decompose and solve shell(s) if there are any: */
  if(nz>=3){
    for(z=1; z<=nz-2; z++){
      rlow = gsl_vector_get(boundlm_z, z-1);
      rhigh = gsl_vector_get(boundlm_z, z);
      alpha = 0.5*(rhigh - rlow);
      beta = 0.5*(rhigh + rlow);
	  for(i=0; i<nn; i++)
		gsl_vector_set(s_shell, i, gsl_matrix_get(sourcelm_zn, z, i));
	  solve_shell(L, alpha, beta, s_shell, f_shell); 
      
	  for(i=0; i<nn; i++){
		gsl_matrix_set(particularlm_zn, z, i, gsl_vector_get(f_shell, i)); 
      }
    }
  }  
  

  /* require continuity of function and derivative: */ 
  solve_continuity(L, boundlm_z, particularlm_zn, homo_coeff);  

  printf("L=%d\n", L);
  print_vector(boundlm_z);
  print_matrix(particularlm_zn);
  print_vector(homo_coeff);

  /* add homogeneous part of solution for kernel */
  alpha = gsl_vector_get(boundlm_z, 0);
  a = gsl_vector_get(homo_coeff, 0);
  if(!(L%2)) /* L is even: */
	chebfit_homo_kernel_even(homo_v, L, a, alpha);
  else /* L is odd: */
	chebfit_homo_kernel_odd(homo_v, L, a, alpha);  
  for(i=0; i<nn; i++)
    gsl_matrix_set(fieldlm_zn, 0, i, gsl_matrix_get(particularlm_zn, 0, i)+gsl_vector_get(homo_v, i)); 
  
  
  /* add homogeneous part of solution for external domain: */
  rlow = gsl_vector_get(boundlm_z, nz-2);
  alpha = -0.5/rlow;
  b = gsl_vector_get(homo_coeff, 2*nz-3);
  chebfit_homo_ext(homo_v, L, b, alpha);  
  for(i=0; i<nn; i++)
	gsl_matrix_set(fieldlm_zn, nz-1, i, gsl_matrix_get(particularlm_zn, nz-1, i)+gsl_vector_get(homo_v, i)); 
  
  
  /* add homogeneous part of solution for shell(s) if there are any: */
  if(nz>=3){
    for(z=1; z<=nz-2; z++){
      rlow = gsl_vector_get(boundlm_z, z-1);
      rhigh = gsl_vector_get(boundlm_z, z);
      alpha = 0.5*(rhigh - rlow);
      beta = 0.5*(rhigh + rlow);
	  a = gsl_vector_get(homo_coeff, 2*z-1);
	  b = gsl_vector_get(homo_coeff, 2*z);
	  chebfit_homo_shell(homo_v, L, a, b, alpha, beta);
	  for(i=0; i<nn; i++)
		gsl_matrix_set(fieldlm_zn, z, i, gsl_matrix_get(particularlm_zn, z, i)+gsl_vector_get(homo_v, i));
    }
  }  
  
  /* free memory */
  gsl_vector_free(s_kernel);
  gsl_vector_free(s_ext);
  gsl_vector_free(s_shell);
  gsl_vector_free(f_kernel);
  gsl_vector_free(f_ext);
  gsl_vector_free(f_shell);
  gsl_vector_free(homo_coeff);
  gsl_vector_free(homo_v);
  gsl_matrix_free(particularlm_zn);
}


void solve_poisson(coeff *f_znlm, coeff *s_znlm, bound_coeff *b_zlm)
{
  int z;
  int n;
  int l;
  int m;
  int nz;
  int nn;
  int nl;
  gsl_vector *boundlm_z;
  gsl_matrix *sourcelm_zn;
  gsl_matrix *fieldlm_zn;
  
  nz = s_znlm->nz;
  nn = s_znlm->nn;
  nl = s_znlm->nl;
  
  boundlm_z = gsl_vector_alloc(nz-1);
  sourcelm_zn = gsl_matrix_alloc(nz, nn);
  fieldlm_zn = gsl_matrix_alloc(nz, nn);
  
  /* solve for real part */
  for(l=0; l<nl; l++){
	for(m=0; m<=l; m++){
	  
	  /* set boundary vector, source matrix for fixed l, m */
	  for(z=0; z<nz-1; z++){
		gsl_vector_set(boundlm_z, z, bound_coeff_get(b_zlm, z, l, m, 0));
		for(n=0; n<nn; n++)
		  gsl_matrix_set(sourcelm_zn, z, n, coeff_get(s_znlm, z, n, l, m, 0));
	  }
	  
	  solve_radial(l, fieldlm_zn, sourcelm_zn, boundlm_z);
	  
	  for(z=0; z<nz-1; z++)
		for(n=0; n<nn; n++)
		  coeff_set(f_znlm, z, n, l, m, 0, gsl_matrix_get(fieldlm_zn, z, n));
	  
	}
  }

  /* solve for immaginary part */
  for(l=0; l<nl; l++){
	for(m=1; m<=l; m++){
	  
	  /* set boundary vector, source matrix for fixed l, m */
	  for(z=0; z<nz-1; z++){
		gsl_vector_set(boundlm_z, z, bound_coeff_get(b_zlm, z, l, m, 1));
		for(n=0; n<nn; n++)
		  gsl_matrix_set(sourcelm_zn, z, n, coeff_get(s_znlm, z, n, l, m, 1));
	  }
	  
	  solve_radial(l, fieldlm_zn, sourcelm_zn, boundlm_z);
	  
	  for(z=0; z<nz-1; z++)
		for(n=0; n<nn; n++)
		  coeff_set(f_znlm, z, n, l, m, 1, gsl_matrix_get(fieldlm_zn, z, n));
	  
	}
  }
  
}



void chebfit_homo_kernel_even(gsl_vector *homo_v, int L, double a, double alpha)
{
  int nn;        /* number of terms in Chebyshev expansion */
  int n;         /* nth coefficient of Chebyshev expansion */
  int i;         /* ith term in summation for coefficient c_n */
  double x_i;    /* Gauss-Lobatto quadrature point */
  double weight;
  double sum;
  
  nn = homo_v->size;
  
  /* Find each of nn coefficients homo_v_n. */
  for(n=0; n<nn; n++){
    sum = 0.0;
    for(i=0; i<nn; i++){
	  /* Evaluate ith term in summation for nth coefficient. */
	  weight = 2.0*cos(PI*n*(nn-1.0-i)/(nn-1))/
		((1.0+delta(0, i)+delta(nn-1, i))*(1.0+delta(0, n)+delta(nn-1, n))*(nn-1));
      x_i = sin(PI*i/(2*(nn-1)));
      sum += a*pow(alpha*x_i, L)*weight;
	}
    gsl_vector_set(homo_v, n, sum);
  }
}


void chebfit_homo_kernel_odd(gsl_vector *homo_v, int L, double a, double alpha)
{
  int nn;        /* number of terms in Chebyshev expansion */
  int n;         /* nth coefficient of Chebyshev expansion */
  int i;         /* ith term in summation for coefficient c_n */
  double x_i;    /* Gauss-Lobatto quadrature point */
  double weight;
  double sum;
  
  nn = homo_v->size;
  
  /* Find each of nn coefficients homo_v_n. */
  for(n=0; n<nn; n++){
    sum = 0.0;
    for(i=0; i<nn; i++){
	  /* Evaluate ith term in summation for nth coefficient. */
	  weight = 2.0*cos(PI*(2*n+1)*(nn-1.0-i)/(2.0*(nn-1)))/
		((1.0+delta(0, i)+delta(nn-1, i))*(1.0+delta(nn-1, n))*(nn-1));
      x_i = sin(PI*i/(2.0*(nn-1)));
      sum += a*pow(alpha*x_i, L)*weight;
	}
    gsl_vector_set(homo_v, n, sum);
  }
}


void chebfit_homo_shell(gsl_vector *homo_v, int L, double a, double b, double alpha, double beta)
{
  int nn;        /* number of terms in Chebyshev expansion */
  int n;         /* nth coefficient of Chebyshev expansion */
  int i;         /* ith term in summation for coefficient c_n */
  double x_i;    /* Gauss-Lobatto quadrature point */
  double weight;
  double sum;
  
  nn = homo_v->size;
  
  /* Find each of nn coefficients homo_v_n. */
  for(n=0; n<nn; n++){
    sum = 0.0;
    for(i=0; i<nn; i++){
	  /* Evaluate ith term in summation for nth coefficient. */
	  weight = 2.0*cos(PI*n*(nn-1.0-i)/(nn-1))/
		((1.0+delta(0, i)+delta(nn-1, i))*(1.0+delta(0, n)+delta(nn-1, n))*(nn-1));
      x_i = -cos(PI*i/(nn-1));
      sum += (a*pow(alpha*x_i+beta, L) + b*pow(alpha*x_i+beta, -L-1))*weight;
	}
    gsl_vector_set(homo_v, n, sum);
  }
}


void chebfit_homo_ext(gsl_vector *homo_v, int L, double b, double alpha)
{
  int nn;        /* number of terms in Chebyshev expansion */
  int n;         /* nth coefficient of Chebyshev expansion */
  int i;         /* ith term in summation for coefficient c_n */
  double x_i;    /* Gauss-Lobatto quadrature point */
  double weight;
  double sum;
  
  nn = homo_v->size;
  
  /* Find each of nn coefficients homo_v_n. */
  for(n=0; n<nn; n++){
    sum = 0.0;
    for(i=0; i<nn; i++){
	  /* Evaluate ith term in summation for nth coefficient. */
	  weight = 2.0*cos(PI*n*(nn-1.0-i)/(nn-1))/
		((1.0+delta(0, i)+delta(nn-1, i))*(1.0+delta(0, n)+delta(nn-1, n))*(nn-1));
      x_i = -cos(PI*i/(nn-1));
      sum += b*pow(alpha*(x_i-1), L+1)*weight;
	}
    gsl_vector_set(homo_v, n, sum);
  }
}
