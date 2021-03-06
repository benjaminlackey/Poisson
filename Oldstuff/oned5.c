/*To compile type: gcc -lm -lgsl -lgslcblas -Wall -pedantic -ansi oned5.c */
/*To write to file type: ./a.out > phigamma2.txt */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define Ne 8 /*Number of terms in matrices for even functions.*/
#define PI 3.141592653589793
#define delta(i, j) (i==j ? 1 : 0)
#define neg1toi(i) (i%2 ? -1 : 1)
#define cosipiby2(i) (i%2 ? 0 : (i%4==0 ? 1 : -1))

void printvector(gsl_vector *v);
void printmatrix(gsl_matrix *m);
void dchebij(gsl_matrix *m);
void inverse(gsl_matrix *m);
void xmatrix(gsl_matrix *m);
double source(double x);
double source_out(double xi);
void chebfit(double a, double b, gsl_vector *c, double (*func)(double));
double chebeval(double a, double b, gsl_vector *c, double x);

int
main (void)
{
  /****************************/
  /*   Declarations           */
  /****************************/
  int i, j;
  double start=-1.0, end=1.0;
  double phisurf; /*potential at boundary of star*/
  double alpha = PI; /*radius of star*/
  double alpha_out = -0.5/alpha; /*external scaling*/ 
  double x;
  int Nall=2*Ne-1;
  gsl_matrix *d = gsl_matrix_alloc (Nall, Nall); /*derivative in T_n basis*/
  gsl_matrix *dsq = gsl_matrix_alloc (Nall, Nall); /*2nd derivative operator*/
  gsl_matrix *dsq_even = gsl_matrix_alloc (Ne, Ne);
  gsl_matrix *inv = gsl_matrix_alloc (Nall, Nall); /*1/x operator*/
  /*gsl_matrix *xop = gsl_matrix_alloc (Nall, Nall);*/ /*x operator*/
  gsl_matrix *L = gsl_matrix_alloc (Nall, Nall); /*linear operator*/ 
  gsl_matrix *L_even = gsl_matrix_alloc (Ne, Ne); /*even part of L*/
  gsl_matrix *invD_even = gsl_matrix_alloc (Ne, Ne);
  gsl_matrix *T = gsl_matrix_alloc (Ne, Ne); /*T_ki = T_i(x_k)*/
  gsl_matrix *TL = gsl_matrix_alloc (Ne, Ne);
  gsl_vector *s = gsl_vector_alloc (Ne); /*source in T_n basis*/
  gsl_vector *s_out = gsl_vector_alloc (Nall);
  gsl_vector *u = gsl_vector_alloc (Ne); /*solution in T_n basis*/
  gsl_vector *u_even_all = gsl_vector_alloc (Nall);
  gsl_vector *u_out = gsl_vector_alloc (Nall); /*solution in T_n basis*/
  int k; /*Used in LU decomposition.  I don't know why.*/
  gsl_permutation * p = gsl_permutation_alloc (Ne); /*Used in LU decomposition.*/
  gsl_permutation * p_out = gsl_permutation_alloc (Nall);
  
  /*gsl_matrix *xsq = gsl_matrix_alloc (Nall, Nall);*/ /**/
	/*gsl_matrix *xfourth = gsl_matrix_alloc (Nall, Nall);*/ /**/
  gsl_matrix *L_out = gsl_matrix_alloc (Nall, Nall); /**/
  gsl_matrix *T_out = gsl_matrix_alloc (Nall, Nall); /**/
  gsl_matrix *TL_out = gsl_matrix_alloc (Nall, Nall); /**/
  
  /*********************************/
  /*      The internal domain.     */
  /*********************************/
  
  /*Make linear differential operator L.*/
  /*L = d^2u/dr^2 + (2/r)du/dr*/
 
  dchebij(d);
  /*printmatrix(d);*/
  inverse(inv);
  /*printmatrix(inv);*/
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, d, d,
		 0.0, dsq);
  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (dsq_even, i, j, gsl_matrix_get(dsq, 2*i, 2*j));
  /*printf("d^2/dx^2 even part:\n");
    printmatrix(dsq_even);*/
  
  /*printf("(1/x)d/dx:\n");*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, inv, d,
		 0.0, L);
  /*printmatrix(L);*/
  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (invD_even, i, j, gsl_matrix_get(L, 2*i, 2*j));
  /* printmatrix(invD_even);*/
  
  
  
  gsl_matrix_scale(L, 2.0);
  gsl_matrix_add(L, dsq); /*WARNING: The order matters here!*/
  /* gsl_matrix_scale(L, 1/(PI*PI));*/
  /*printmatrix(L);*/
  
  /*xmatrix(xop, Nall);*/
  /*printmatrix(xop);*/

  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (L_even, i, j, gsl_matrix_get(L, 2*i, 2*j));
  /*printmatrix(L_even);*/
  
  /*Make colocation matrix T_ki = T_i(x_k) where x_k = sin(k*PI/(2*ORDER)).*/
  for(i=0; i<Ne; i++)
    for(j=0; j<Ne; j++)
      gsl_matrix_set (T, i, j, neg1toi(j)*cos(i*j*PI/(Ne-1)));
  /*printmatrix(T);*/
  
  /*Make operator T.L*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, T, L_even,
		 0.0, TL);
  /*printmatrix(TL);*/
  
  /*Make source vector.*/
  for(i=0; i<Ne; i++){
    gsl_vector_set(s, i, alpha*alpha*source(sin(i*PI/(2*(Ne-1)))));
    /*printf("%f\n", sin(i*PI/(2*(Ne-1))));*/
  }
  /*printvector(s);*/
  
  /*Set boundary conditions in TL and s.*/
  for(i=0; i<Ne; i++){
    gsl_matrix_set (TL, 0, i, neg1toi(i));
  }
  gsl_vector_set(s, 0, -2);
  /*printmatrix(TL);
    printvector(s);*/
  
  /*Solve for u~ in T.L.u~ = s~.*/
  gsl_linalg_LU_decomp (TL, p, &k); 
  gsl_linalg_LU_solve (TL, p, s, u); 
  /*printvector(u);*/
  

  gsl_vector_set_zero(u_even_all);
  for(i=0; i<Ne; i++){
	gsl_vector_set(u_even_all, 2*i, gsl_vector_get(u, i));
  }
  
  /*for(x=0.0; x<=end; x+=0.01)
    printf("%f\t%f\n", x, chebeval(start, end, u_even_all, x));*/
  for(x=0.0; x<=end; x+=0.01)
	printf("%f\t%.18e\n", x, chebeval(start, end, u_even_all, x)-(-1-source(x)/(PI*PI)));
  
  /*********************************/
  /*      The external domain.     */
  /*********************************/
  
/*   xmatrix(xop, Nall); */
/*   printmatrix(xop, Nall, Nall); */

/*   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,  */
/* 		 1.0, xop, xop, */
/* 		 0.0, xsq); */
/*   printmatrix(xsq, Nall, Nall); */
/*   /\*x.x.x.x :*\/ */
/*   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,  */
/* 		 1.0, xsq, xsq, */
/* 		 0.0, xfourth); */
/*   printmatrix(xfourth, Nall, Nall); */
/*   /\*Now L = x^4(d^2/dx^2) :*\/ */
/*   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,  */
/* 		 1.0, xfourth, dsq, */
/* 		 0.0, L_out); */
/*   printmatrix(L_out, Nall, Nall); */

  gsl_matrix_memcpy(L_out, dsq); /*Set L_out to dsq.*/
  /*printmatrix(L_out);*/
  
  /*Make colocation matrix T_out.*/
  for(i=0; i<Nall; i++)
    for(j=0; j<Nall; j++)
      gsl_matrix_set (T_out, i, j, neg1toi(j)*cos(i*j*PI/(Nall-1)));
  /*printmatrix(T_out);*/
  
  /*Make operator T.L*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 
		 1.0, T_out, L_out,
		 0.0, TL_out);
  /*printmatrix(TL_out);*/
  
  /*Make source vector.*/
  for(i=0; i<Nall; i++){
    gsl_vector_set(s_out, i, source_out(-cos(i*PI/(Nall-1)))/(alpha_out*alpha_out));
    /*printf("%f\n", -cos(i*PI/(Nall-1)));*/
  }
  /*printvector(s_out);*/
  
  /*Set boundary conditions in TL and s.*/
  for(i=0; i<Nall; i++){
    gsl_matrix_set (TL_out, 0, i, neg1toi(i));
    gsl_matrix_set (TL_out, Nall-1, i, 1);
  }
  gsl_vector_set(s_out, 0, -1);
  gsl_vector_set(s_out, Nall-1, 0);
  
  /*Solve for u~ in T.L.u~ = s~.*/
  gsl_linalg_LU_decomp (TL_out, p_out, &k); 
  gsl_linalg_LU_solve (TL_out, p_out, s_out, u_out); 
  /*printvector(u_out);*/
  
  /*for(x=-1.0; x<=end; x+=0.01)
    printf("%f\t%f\n", x, chebeval(start, end, u_out, x));*/
  for(x=-1.0; x<=end; x+=0.01)
    printf("%f\t%.18e\n", x, chebeval(start, end, u_out, x)-0.5*(x-1));
  
  /*********************************/
  /*      Free some memory.        */
  /*********************************/
  
  gsl_matrix_free(d);
  gsl_matrix_free(dsq);
  gsl_matrix_free(dsq_even);
  gsl_matrix_free(inv);
  gsl_matrix_free(L);
  gsl_matrix_free(L_even);
  gsl_matrix_free(invD_even);
  gsl_matrix_free(T);
  gsl_matrix_free(TL);
  gsl_vector_free(s);
  gsl_vector_free(s_out);
  gsl_vector_free(u);
  gsl_vector_free(u_even_all);
  gsl_vector_free(u_out);
  gsl_matrix_free(L_out);
  gsl_matrix_free(T_out);
  gsl_matrix_free(TL_out);

  return 0;
}

/*************************************/
/*    Matrix operator functions.     */
/*************************************/

/*Returns [0...(n-1)]x[0...(n-1)] matrix for derivative of Chebyshev basis.*/  
void dchebij(gsl_matrix *m)
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

/*RECHECK THIS MATRIX.*/
void xmatrix(gsl_matrix *m)
{
  int n; /*nXn matrix*/
  int i; /*row*/ 
  
  n=(*m).size1;

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

void inverse(gsl_matrix *m)
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
  for(i=1; i<n; i++){
    sign=1;
    for(j=i+1; j<n; j+=2){
      gsl_matrix_set (m, i, j, sign*2);
      sign*=-1;
    }
  }
}


/*******************************/
/*      Source functions.      */
/*******************************/

double source(double xi)
{
  return 1;

  /*if(xi<0.01)
    return 1-pow(xi*PI, 2)/6.0+pow(xi*PI, 4)/120.0;
  else
    return sin(PI*xi)/(PI*xi);*/
}

/*rho(xi)/(xi-1)^4*/
double source_out(double xi)
{
  return 
}


/*******************************************************/
/* Going between functions and Chebyshev coefficients. */
/*******************************************************/

/*Takes a function and finds N Chebyshev coefficients      */
/*upto order N-1 using Chebyshev-Gauss-Lobatto quadrature. */
void chebfit(double a, double b, gsl_vector *c, double (*func)(double))
{
  int N;
  int k;            /*kth coefficient of Chebyshev expansion*/
  int j;            /*jth term in summation for coefficient c_k*/
  double middle;     /*midpoint of range [a, b]*/ 
  double half_width; /*half width of [a, b]*/
  double x_j;        /*Gauss-Lobatto quadrature point*/
  double sum;       /*sum of terms in summation for c_k*/
  double fac;        /*overall factor in front of summation for c_k*/
  
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

/*Takes N coefficients of Chebyshev polynomials upto order N-1*/
/*and uses Clenshaw's recurrance formula to determine the     */
/*function value at x.                                        */
double chebeval(double a, double b, gsl_vector *c, double x)
{ 
  double t; /*change of variable*/ 
  int N;
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

void printmatrix(gsl_matrix *m)
{
  int rows, columns;
  int i, j;

  rows = (*m).size1;
  columns = (*m).size2;
  
  for (i = 0; i < rows; i++){
    for (j = 0; j < columns; j++)
      printf ("%4g  ", gsl_matrix_get (m, i, j));
    printf ("\n");
  } 
  printf("\n");
}

void printvector(gsl_vector *v)
{
  int n;
  int i;

  n = (*v).size;
  
  for (i = 0; i < n; i++){
    printf ("%4g  ", gsl_vector_get (v, i));
  } 
  printf("\n\n");
}
