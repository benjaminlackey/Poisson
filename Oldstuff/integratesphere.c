/* To compile type: gcc -g -lm -lgsl -lgslcblas -Wall -pedantic -ansi integratesphere.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

typedef struct 
{
  int nz;
  int nn;
  int nl;
  double *data;
} coeff;

typedef struct
{
  int nz;
  int nl;
  double *data;
} bound_coeff;

/* constants and macros */
#define PI 3.141592653589793
#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
#define cosipiby2(i) ((i)%2 ? 0 : ((i)%4==0 ? 1 : -1)) /* \cos(i\pi/2) */

/* function prototypes */
/* tested: */
void legendre_zero_weight(gsl_vector *zeros, gsl_vector *weights);

coeff *coeff_alloc(int z, int n, int L);
void coeff_set(coeff *c, int z, int n, int L, int m, int i, double x);
double coeff_get(coeff *c, int z, int n, int L, int m, int i);
void print_coeff(coeff *c);

bound_coeff *bound_coeff_alloc(int z, int L);
void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x);
double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i);
void print_bound_coeff(bound_coeff *b);
void print_vector(gsl_vector *v);
/* not tested: */

double source(double r, double theta, double phi);
double boundary(int zone, double theta, double phi);
double p_lm(int l, int m, double x);
double ytheta_lm(int l, int m, double x);
void decompose_boundary(bound_coeff *b_zlm, 
						int nphi,
						int nx,
						double (*boundary)(int, double, double));
/*void decompose_source(coeff *c_znlm, 
					  bound *b_zlm, 
					  double (*source)(double, double, double)); */


int main(void)
{
  bound_coeff *b_zlm = bound_coeff_alloc(2, 7);
  gsl_vector *zeros = gsl_vector_alloc(16);
  gsl_vector *weights = gsl_vector_alloc(16);
  int i;
  double sum;
  
  legendre_zero_weight(zeros, weights);
  print_vector(zeros);
  print_vector(weights);
  sum=0;
  for(i=0; i<16; i++)
	sum += gsl_vector_get(weights, i);
  printf("%f\n\n", sum);

  decompose_boundary(b_zlm, 100, 100, boundary);
  
  print_bound_coeff(b_zlm);

  printf("%.18e\n", ytheta_lm(0, 0, cos(3.1)));

  return 0;
}


double boundary(int zone, double theta, double phi)
{
  /*if(zone==0)
	return sqrt(1/(4*PI)) + sqrt(3/(4*PI))*cos(theta) - sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi); 
  else
  return 7*sqrt(1/(4*PI)) + 3*sqrt(3/(4*PI))*cos(theta) - 2*sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi);*/
  /*return 1;*/
  if(zone==0)
	return sin(theta*phi); 
  else
	return 7*sqrt(1/(4*PI)) + 3*sqrt(3/(4*PI))*cos(theta) - 2*sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi);
}


double source(double r, double theta, double phi)
{
  return 1; 
}

/*************************************************************/
/* For each boundary decompose it into spherical harmonics.  */
/* The boundary is a real function so negative m is ignored. */
/*************************************************************/
void decompose_boundary(bound_coeff *b_zlm, 
						int nphi,
						int nx,
						double (*boundary)(int, double, double))
{
  int z; /* current boundary */
  int nz; /* number of boundaries */
  int nl; /* number of multipole moments lmax=nl-1 */ 
  int l; /* current multipole moment */
  int m; /* current azimuthal part */
  int i; /* iterator for phi integration */
  int j; /* iterator for theta integration */
  double phisum;
  double thetasum;
  double x_j; /* Legendre polynomial zero */
  double w_j; /* weight */
  gsl_vector *zeros = gsl_vector_alloc(nx);
  gsl_vector *weights = gsl_vector_alloc(nx);
  nz = b_zlm->nz;
  nl = b_zlm->nl;

  /*printf("%d\t%d\n", nz, nl);*/
  /* find coefficients for each boundary */
  for(z=0; z<nz; z++){
	for(l=0; l<nl; l++){
	  /* integrate real part */
	  for(m=0; m<=l; m++){
		/* sum up terms for theta part */
		thetasum = 0;
		legendre_zero_weight(zeros, weights);
		for(j=0; j<nx; j++){
		  x_j = gsl_vector_get(zeros, j);
		  w_j = gsl_vector_get(weights, j);
		  /* sum up terms for phi part */
		  phisum = 0;
		  for(i=0; i<nphi; i++){
			phisum += cos(2*PI*m*i/nphi)*boundary(z, acos(x_j), 2*PI*i/nphi);
		  }
		  thetasum += w_j*ytheta_lm(l, m, x_j)*phisum; 
		}
		thetasum *= 2*PI/nphi;
		/* thetasum is now coefficient of boundary for given z, l, m */
		bound_coeff_set(b_zlm, z, l, m, 0, thetasum);
	  }
	  
	  /* integrate imaginary part */
	  for(m=1; m<=l; m++){
		/* sum up terms for theta part */
		thetasum = 0;
		legendre_zero_weight(zeros, weights);
		for(j=0; j<nx; j++){
		  x_j = gsl_vector_get(zeros, j);
		  w_j = gsl_vector_get(weights, j);
		  /* sum up terms for phi part */
		  phisum = 0;
		  for(i=0; i<nphi; i++){
			phisum += -sin(2*PI*m*i/nphi)*boundary(z, acos(x_j), 2*PI*i/nphi);
		  }
		  thetasum += w_j*ytheta_lm(l, m, x_j)*phisum;
		}
		thetasum *= 2*PI/nphi;
		/* thetasum is now coefficient of boundary for given z, l, m */
		bound_coeff_set(b_zlm, z, l, m, 1, thetasum);
	  }
	  
	}	
  }
}


/****************************************************************************/
/* Decompose a function into spherical harmonics and chebyshev polynomials. */
/****************************************************************************/
void decompose_boundary(bound_coeff *b_zlm, 
						int nphi,
						int nx,
						double (*boundary)(int, double, double))
{




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
	

/***********************************************/
/* Calculate the associated Legendre functions */
/***********************************************/
double p_lm(int l, int m, double x)
{
  double pmm=1.0; /*p00=1*/
  double a;
  double fact=1.0;
  int i;
  double pmmp1;
  int ll;
  double pmmp2;

  if(m>0) { /*find pmm for m>0*/
    a= -sqrt((1.0-x)*(1.0+x));
    for(i=1; i<=m; i++) {
      pmm*= fact*a;
      fact+=2.0;
    }
  }
  
  if(l==m) {
    return pmm;
  } else { /*find plm for l>m*/
	pmmp1=x*(2*m+1)*pmm;
    if(l==m+1) {
      return pmmp1; /*here l=m+1*/
    } else { /*here l>=m+2*/
      for(ll=m+2; ll<=l; ll++) {
		pmmp2=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m); /*using recursion relation*/
		pmm=pmmp1;
		pmmp1=pmmp2;
      }
      return pmmp2;
    }
  }
}

/*******************************************************/
/* Calculate the theta part of the spherical harmonic  */
/* Y_{lm} without the exp(i*m*phi)) factor.            */
/*******************************************************/
double ytheta_lm(int l, int m, double x)
{
  int i;
  double fact;
  double ylmabs;

  if(l==0 && m==0) {
    return 1/sqrt(4*PI);
  } else {
    fact=1.0;
    for(i=l-m+1; i<=l+m; i++)
      fact*= i;
    
    ylmabs=sqrt((2*l+1)/(4*PI*fact))*p_lm(l, m, x);

    if(m<0) {
      if((-m)%2) {
	return -ylmabs;
      } else {
	return ylmabs;
      }
    } else {
      return ylmabs;
    }
  }
}

/*-------------------------------------------------------------------*/


/*****************************************/
/* Allocate memory to store coefficients */
/* and return a pointer to the memory.   */
/*****************************************/
coeff *coeff_alloc(int nz, int nn, int nl)
{
  /* z goes from 0 to nz-1 */
  /* n goes from 0 to nn-1 */
  /* l goes from 0 to lmax=nl-1 */
  /* m goes from -l to l for l */
  
  coeff *c;
  double *data;
  
  /* sizeof(coeff) is the number of bytes needed to store 3 integers
     and the address of data.  Memory for the actual data is not
     allocated here: */
  c = (coeff *)malloc(sizeof(coeff));
  
  /* Memory for the data is allocated here: */
  /* If -l <= m <= l for each l, there are (lmax+1)^2 elements
	 for each n and z (or nl^2 elements where nl=lmax+1). */
  data = (double *)malloc(sizeof(double)*nz*nn*nl*nl);
  
  c->nz = nz;
  c->nn = nn;
  c->nl = nl;
  c->data = data;
  
  return c;
}


/**************************************************/
/* Allocate memory to store boundary coefficients */
/* and return a pointer to the memory.            */
/**************************************************/
bound_coeff *bound_coeff_alloc(int nz, int nl)
{
  /* z goes from 0 to nz-1 */
  /* l goes from 0 to lmax=nl-1 */
  /* m goes from -l to l for l */
  
  bound_coeff *b;
  double *data;
  
  b = (bound_coeff *)malloc(sizeof(bound_coeff));
  
  data = (double *)malloc(sizeof(double)*nz*nl*nl);
  
  b->nz = nz;
  b->nl = nl;
  b->data = data;
  
  return b;
}


/********************************/
/* Set coefficient c_znlm to x. */
/********************************/
/*extern inline */
void coeff_set(coeff *c, int z, int n, int L, int m, int i, double x)
{ 
  /* do error checking */
  if((z >= c->nz)||(n >= c->nn)||(L >= c->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in coeff_set\n");
  
  /* index = z*nn*nl^2 + n*nl^2 + L^2 + 2m + i + (delta_{0m} - 1) */
  /*       = ((z*nn + n)*nl^2 + L^2 + 2m + i + (delta_{0m} - 1) */
  c->data[(z*c->nn + n)*c->nl*c->nl + L*L + 2*m + i + delta(0, m) - 1] = x ;
}


/********************************************/
/* Set coefficient b_zlm for boundary to x. */
/********************************************/
/*extern inline */
void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x)
{ 
  /* do error checking */
  if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_set\n");
  
  /* index = z*nl^2 + n*nl^2 + L^2 + 2m + i+ (delta_0m - 1) */
  b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1] = x ;
}


/*****************************/
/* Get value of c_znlm.      */
/*****************************/
/*extern inline */
double coeff_get(coeff *c, int z, int n, int L, int m, int i)
{
  /* do error checking */
  if((z >= c->nz)||(n >= c->nn)||(L >= c->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in coeff_get\n");
  
  return c->data[(z*c->nn + n)*c->nl*c->nl + L*L + 2*m + i + delta(0, m) - 1];
} 

/*****************************/
/* Get value of b_zlm.      */
/*****************************/
/*extern inline */
double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i)
{
  /* do error checking */
  if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_get\n");
  
  return b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1];
}


/*******************************/
/* Print all the coefficients. */
/*******************************/
void print_coeff(coeff *c){
  int nz = c->nz;
  int nn = c->nn;
  int nl = c->nl;
  int z, n, L, m;
  
  for(z=0; z<nz; z++) {
    printf("zone %d:\n", z);
    for(n=0; n<nn; n++) {
      printf("radial basis %d:\n\t", n);
	  printf("m=0\t");
      for(m=1; m<nl; m++)
		printf("m=%d\t\t", m);
      printf("\n");
	  printf("L=0:\t");
	  printf("(%f, -)\n", coeff_get(c, z, n, 0, 0, 0));
      for(L=1; L<nl; L++) {
		printf("L=%d:\t", L);
		printf("(%f, -)\t", coeff_get(c, z, n, L, 0, 0));
		for(m=1; m<=L; m++){
		  printf("(%f, %f)\t", coeff_get(c, z, n, L, m, 0), coeff_get(c, z, n, L, m, 1));
		}
		printf("\n");
      }
    }
  }
  printf("\n");
}


/****************************************/
/* Print all the boundary coefficients. */
/****************************************/
void print_bound_coeff(bound_coeff *b){
  int nz = b->nz;
  int nl = b->nl;
  int z, L, m;
  
  for(z=0; z<nz; z++) {
    printf("zone %d:\n", z);
	printf("\tm=0\t");
	for(m=1; m<nl; m++)
	  printf("m=%d\t\t", m);
	printf("\n");
	printf("L=0:\t");
	printf("(%f, -)\n", bound_coeff_get(b, z, 0, 0, 0));
	for(L=1; L<nl; L++) {
	  printf("L=%d:\t", L);
	  printf("(%f, -)\t", bound_coeff_get(b, z, L, 0, 0));
	  for(m=1; m<=L; m++){
		printf("(%f, %f)\t", bound_coeff_get(b, z, L, m, 0), bound_coeff_get(b, z, L, m, 1));
	  }
	  printf("\n");
	}
  }
  printf("\n");
}

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
