/* gcc -g -lm -lgsl -lgslcblas -Wall -pedantic -ansi chebyshevtest.c print.c coefficients.c poisson.h */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include "poisson.h"

void chebfit(gsl_vector *c, double (*func)(double));
void chebfit_even(gsl_vector *c, double (*func)(double));
void chebfit_odd(gsl_vector *c, double (*func)(double));

double function(double x);
double evenfunction(double x);
double oddfunction(double x);

int main(void)
{
  gsl_vector *c = gsl_vector_alloc(16);
  gsl_vector *homo_v = gsl_vector_alloc(8);

  chebfit(c, function);

  print_vector(c);
  
  chebfit_even(c, evenfunction);

  print_vector(c);

  chebfit_odd(c, oddfunction);

  print_vector(c);
  
  chebfit_homo_kernel_even(homo_v, 2, 2, 2);
  print_vector(homo_v);
  
  chebfit_homo_kernel_odd(homo_v, 3, 2, 2);
  print_vector(homo_v);
  return 0;
}

double function(double x)
{
  return 1 - 2*x - 6*x*x + 4*x*x*x + 8*x*x*x*x;
}

double evenfunction(double x)
{
  return 1 - 20*x*x + 120*pow(x, 4) - 224*pow(x, 6) + 128*pow(x, 8);
  /*return 3 - 54*x*x + 488*pow(x, 4) - 1056*pow(x, 6) + 640*pow(x, 8);*/
}

double oddfunction(double x)
{
  /*return 5*x - 80*pow(x, 3) + 336*pow(x, 5) - 512*pow(x, 7) + 256*pow(x, 9);*/
  return 95*x - 984*pow(x, 3) + 2960*pow(x, 5) - 3584*pow(x, 7) + 1536*pow(x, 9);
}


/*************************************************************/
/* Takes a function and finds N Chebyshev coefficients       */
/* upto order N-1 using Chebyshev-Gauss-Lobatto quadrature.  */
/*************************************************************/
void chebfit(gsl_vector *c, double (*func)(double))
{
  int nn;        /* number of terms in Chebyshev expansion */
  int n;         /* nth coefficient of Chebyshev expansion */
  int i;         /* ith term in summation for coefficient c_n */
  double x_i;    /* Gauss-Lobatto quadrature point */
  double weight;
  double sum;
  
  nn = c->size;
  
  /* Find each of nn coefficients c_n. */
  for(n=0; n<nn; n++){
    sum = 0.0;
    for(i=0; i<nn; i++){
	  /* Evaluate ith term in summation for nth coefficient. */
	  weight = 2.0*cos(PI*n*(nn-1.0-i)/(nn-1))/
		((1.0+delta(0, i)+delta(nn-1, i))*(1.0+delta(0, n)+delta(nn-1, n))*(nn-1));
      x_i = -cos(PI*i/(nn-1));
      sum += func(x_i)*weight;
	}
    gsl_vector_set(c, n, sum);
  }
}


void chebfit_even(gsl_vector *c, double (*func)(double))
{
  int nn;        /* number of terms in Chebyshev expansion */
  int n;         /* nth coefficient of Chebyshev expansion */
  int i;         /* ith term in summation for coefficient c_n */
  double x_i;    /* Gauss-Lobatto quadrature point */
  double weight;
  double sum;
  
  nn = c->size;
  
  /* Find each of nn coefficients c_n. */
  for(n=0; n<nn; n++){
    sum = 0.0;
    for(i=0; i<nn; i++){
	  /* Evaluate ith term in summation for nth coefficient. */
	  weight = 2.0*cos(PI*n*(nn-1.0-i)/(nn-1))/
		((1.0+delta(0, i)+delta(nn-1, i))*(1.0+delta(0, n)+delta(nn-1, n))*(nn-1));
      x_i = sin(PI*i/(2*(nn-1)));
      sum += func(x_i)*weight;
	}
    gsl_vector_set(c, n, sum);
  }
}


void chebfit_odd(gsl_vector *c, double (*func)(double))
{
  int nn;        /* number of terms in Chebyshev expansion */
  int n;         /* nth coefficient of Chebyshev expansion */
  int i;         /* ith term in summation for coefficient c_n */
  double x_i;    /* Gauss-Lobatto quadrature point */
  double weight;
  double sum;
  
  nn = c->size;
  
  /* Find each of nn coefficients c_n. */
  for(n=0; n<nn; n++){
    sum = 0.0;
    for(i=0; i<nn; i++){
	  /* Evaluate ith term in summation for nth coefficient. */
	  weight = 2.0*cos(PI*(2*n+1)*(nn-1.0-i)/(2.0*(nn-1)))/
		((1.0+delta(0, i)+delta(nn-1, i))*(1.0+delta(nn-1, n))*(nn-1));
      x_i = sin(PI*i/(2.0*(nn-1)));
      sum += func(x_i)*weight;
	}
    gsl_vector_set(c, n, sum);
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
