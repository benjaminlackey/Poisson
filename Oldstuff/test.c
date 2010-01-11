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
void legendre_zero_weight(gsl_vector *zeros, gsl_vector *weights);
void print_vector(gsl_vector *v);

int main(void)
{
  gsl_vector *zeros = gsl_vector_alloc(16);
  gsl_vector *weights = gsl_vector_alloc(16);
  
  legendre_zero_weight(zeros, weights);
  print_vector(zeros);
  print_vector(weights);

  return 0;
}


/*************************************************************/
/* Find the zeros x_j of the nth Legendre polynomial         */
/* and the corresponding weights w_j for Gaussian quadrature */
/*************************************************************/
void legendre_zero_weight(gsl_vector *zeros, gsl_vector *weights)
{
  int n;              /* order of polynomial and number of zeros */
  double acc=1.0e-15; /* absolute accuracy */
  double pn;          /* nth Legendre polynomial P_n(x) */
  double pnminus1;    /* (n-1)th Legendre polynomial P_{n-1}(x) */ 
  double pnminus2;    /* (n-2)th Legendre polynomial P_{n-2}(x) */ 
  double pnprime;     /* derivative of nth Legendre polynomial P_n'(x) */
  double x;           /* current guess for zero */
  double xold;        /* old guess for zero */
  int i;              /* ith zero to the right of the origin */
  int j; 
  int m;              /* number of zeros to the right of the origin */
  
  n = zeros->size;
  m = (n+1)/2;
  for(i=1; i<=m; i++) {
	/* Initial guess for the zero */
	x = cos(PI*(i-0.25)/(n+0.5));

	/* Start Newton's method to find the zero */
	do {
	  /* Evaluate nth Legendre polynomial at point x.      */
	  /* Use the recursion relation:                       */
	  /* n*P_n(x) = (2n-1)*x*P_{n-1}(x) - (n-1)*P_{n-2}(x) */
	  pn = 1.0;
	  pnminus1 = 0.0;
	  for(j=1; j<=n; j++) {
		pnminus2 = pnminus1;
		pnminus1 = pn;
		pn=((2.0*j-1.0)*x*pnminus1 - (j-1.0)*pnminus2)/j;
	  }
	 
	  /* Evaluate derivative of nth Legendre polynomial at point x. */
	  /* Use the relation:                                          */
	  /* (1-x^2)*P_n'(x) = -n*x*P_n(x) + n*P_{n-1}(x)               */
	  pnprime = n*(-x*pn + pnminus1)/(1.0-x*x);
	  
	  /* Make next guess for the zero */
	  xold = x;
	  x = xold - pn/pnprime;
	} while(fabs(x-xold) > acc);

	/* Set vectors for zeros and weights */ 
	gsl_vector_set(zeros, n-i, x);
	gsl_vector_set(zeros, i-1, -x);
	gsl_vector_set(weights, n-i, 2/((1-x*x)*pnprime*pnprime));
	gsl_vector_set(weights, i-1, 2/((1-x*x)*pnprime*pnprime));
  }
}


void print_vector(gsl_vector *v)
{
  int n;
  int i;
  
  n = v->size;
  
  for (i = 0; i < n; i++){
    printf ("%.18e\t", gsl_vector_get (v, i));
  } 
  printf("\n\n");
}
