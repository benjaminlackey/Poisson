/* To compile type: gcc -lm -lfftw3 fouriertoplm.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* fftw header */
# include <fftw3.h>

#define PI 3.141592653589793

#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */

double p_lm(int l, int m, double x);
double ytheta_lm(int l, int m, double x);

int main (void)
{
  int l=3;
  int m=2;
  int i;
  double x_i;
  double *in;
  double *in2;
  int n = 10 /*2*l+1*/;
  int nc;
  fftw_complex *out;
  fftw_plan plan_backward;
  fftw_plan plan_forward;
  
  /* Set up an array to hold the data, and assign the data. */
  in = fftw_malloc ( sizeof ( double ) * n );
  
  for ( i = 0; i < n; i++ )
	{
	  x_i = 2*PI*i/n;
	  /*in[i] = p_lm(l, m, cos(x_i));*/
	  in[i] = 1 + cos(x_i) + sin(x_i) + cos(5*x_i) + sin(4*x_i);
	}
  
  printf ( "\n" );
  printf ( "  Input Data:\n" );
  printf ( "\n" );
  
  for ( i = 0; i < n; i++ )
	{
	  printf ( "  %4d  %12f\n", i, in[i] );
	}
  /*
	Set up an array to hold the transformed data,
	get a "plan", and execute the plan to transform the IN data to
	the OUT FFT coefficients.
  */
  nc = ( n / 2 ) + 1;
  
  out = fftw_malloc ( sizeof ( fftw_complex ) * nc );
  
  plan_forward = fftw_plan_dft_r2c_1d ( n, in, out, FFTW_ESTIMATE );
  
  fftw_execute ( plan_forward );
  
  printf ( "\n" );
  printf ( "  Output FFT Coefficients:\n" );
  printf ( "\n" );
  
  for ( i = 0; i < nc; i++ )
	{
	  printf ( "  %4d  %12f  %12f\n", i, out[i][0]/n, out[i][1]/n );
	}
  /*
	Set up an arrray to hold the backtransformed data IN2,
	get a "plan", and execute the plan to backtransform the OUT
	FFT coefficients to IN2.
  */
  in2 = fftw_malloc ( sizeof ( double ) * n );
  
  plan_backward = fftw_plan_dft_c2r_1d ( n, out, in2, FFTW_ESTIMATE );
  
  fftw_execute ( plan_backward );
  
  printf ( "\n" );
  printf ( "  Recovered input data divided by N:\n" );
  printf ( "\n" );
  
  for ( i = 0; i < n; i++ )
	{
	  printf ( "  %4d  %12f\n", i, in2[i]/n /*( double ) ( n )*/ );
	}
  /*
	Release the memory associated with the plans.
  */
  fftw_destroy_plan ( plan_forward );
  fftw_destroy_plan ( plan_backward );
  
  fftw_free ( in );
  fftw_free ( in2 );
  fftw_free ( out );
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
