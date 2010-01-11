#include <stdio.h>
#include <math.h>

#define PI 3.141592653589793

double ftheta(double x);
double p_n(int n, double x);
double p_nprime(int n, double x);
double p_nroot(int n, double a, double b);
void tabulateroots(int n, double *roottable);
double thetaint(int n, double *roottable);
double rephiint(int m);
double imphiint(int m);
double refphi(double phi);
double imfphi(double phi);
double p_lm(int l, int m, double x);
double ytheta_lm(int l, int m, double x);

main()
{
  int n=10;
  /*make a list of the roots of p_n(n, x)*/
  double list[1000];
  tabulateroots(n, &list[0]); 
  /*print the roots*/
  printf("The roots of the %dth order Legendre polynomial are:\n", n);
  int i;
  for(i=0; i<n; i++)
    printf("%f\n", list[i]);

  /*evaluate the integral y*43 times y43*/

  double re;
  double im;

  re=thetaint(n, &list[0])*rephiint(n);
  im=thetaint(n, &list[0])*imphiint(n);

  printf("\nint(Y*43 Y43) = %16.13f + i%16.13f\n", re, im); 
}


/*********************************************
 * Functions to test the integration scheme. *
 *********************************************/

surf(float theta, float phi)
{
  return 3*pow(cos(theta), 2)-1 + pow(sin(theta), 2)*cos(theta)*cos(2*phi);
}


/*Function for the theta coordinate.  x=cos(theta).*/
double ftheta(double x)
{
  return ytheta_lm(4, 3, x)*ytheta_lm(4, 3, x);

  /*return x*(1-x)*exp(-x);*/
  /*return 7*pow(x, 6)-3*pow(x, 5)+2*pow(x, 4)+5*pow(x, 3)-7*pow(x, 2)+x-2;*/
}

/*Function for the real part of the phi coordinate.*/
double refphi(double phi)
{
  return cos((3-3)*phi);
}

/*Function for the imaginary part of the phi coordinate.*/
double imfphi(double phi)
{
  return sin((3-3)*phi);
}


/*********************************************
 * Gauss-Legendre integration of theta part. *
 *********************************************/

/*Evaluate the nth Legendre polynomial with a recursion relationship.*/
double p_n(int n, double x)
{
  if (n==0)
    return 1;
  else if (n==1) 
    return x;
  else
    return ((2*n-1)*x*p_n(n-1, x)-(n-1)*p_n(n-2, x))/n;
}


/*Evaluate the derivative of the nth Legendre polynomial.*/
double p_nprime(int n, double x)
{
  if(n==0) /*special case where formula below doesn't work*/
    return 0;
  
  return (-n*x*p_n(n, x)+n*p_n(n-1, x))/(1-x*x);
}


/*Find roots with bisection.*/    
double p_nroot(int n, double xl, double xh)
{
  double acc=1.0e-13; /*absolute accuracy*/
  int maxit=50; /*maximum number of iterations*/
  double xm;
  double yxl=p_n(n, xl); 
  double yxh=p_n(n, xh);
  double yxm; 
  
  if(yxl*yxh>0){
    printf("%f and %f do not bracket the root.\n", xl, xh);
    return 0.5*(xl+xh); /*return something*/
  }
  
  int i;
  for(i=1; i<=maxit; i++) {
    xm=0.5*(xl+xh); /*mid point*/ 
    yxm=p_n(n, xm);
    
    if(yxl*yxm<=0) /*Root is bracketed by xl and xm.*/
      xh=xm;
    else /*Root is bracketed by xm and xh.*/
      xl=xm;
    
    if(xh-xl<=acc)
      return xm;
  }
  printf("Number of iterations has exceeded %d.\n", maxit);
  printf("Accuracy is less than %e.\n", acc);
  return xm;
}


/*Find all the roots for a Legendre polynomial of order n.*/
void tabulateroots(int n, double *roottab)
{
  int max;
  max=n/2;
  
  if(n%2)
    roottab[max]=0.0; /*the middle root is 0 for odd n*/
  
  double a, b, root;
  int m;
  for(m=max; m>0; m--) {
    a=cos(2*m*PI/(2*n+1));
    b=cos((2*m-1)*PI/(2*n+1));
    root=p_nroot(n, a, b); 
    roottab[n-m]=root;
    roottab[m-1]= -root;
  }
  
}


/*Integrate a function by summing the roots and the weights.*/
double thetaint(int n, double *roottable)
{
  double sum=0;
  double w;
  int i;
  for(i=0; i<n; i++){
    w=2.0/((1-roottable[i]*roottable[i])
	   *p_nprime(n, roottable[i])*p_nprime(n, roottable[i]));
    sum+= w*ftheta(roottable[i]);
  }
  return sum;
}

/*********************************************
 * Integration of phi part.                  *
 *********************************************/


/*Integrate over the phi direction (real part)*/
double rephiint(int m)
{
  double sum=0;
  double w;
  int i;
  for(i=0; i<m; i++) {
    w= 2*PI/m;
    sum+= w*refphi(2*PI*i/m);
  }
  return sum;
}


/*Integrate over the phi direction (imaginary part)*/
double imphiint(int m)
{
  double sum=0;
  double w;
  int i;
  for(i=0; i<m; i++) {
    w= 2*PI/m;
    sum+= w*imfphi(2*PI*i/m);
  }
  return sum;
}

/*********************************************
 * Calculate special functions.              *
 *********************************************/

/*Calculate the associated Legendre functions*/
double p_lm(int l, int m, double x)
{
  double pmm=1.0; /*p00=1*/
  double a;
  double fact=1.0;

  if(m>0) { /*find pmm for m>0*/
    a= -sqrt((1.0-x)*(1.0+x));
    int i;
    for(i=1; i<=m; i++) {
      pmm*= fact*a;
      fact+=2.0;
    }
  }

  if(l==m) {
    return pmm;
  } else { /*find plm for l>m*/
    double pmmp1=x*(2*m+1)*pmm;
    if(l==m+1) {
      return pmmp1; /*here l=m+1*/
    } else { /*here l>=m+2*/
      int ll;
      double pmmp2;
      for(ll=m+2; ll<=l; ll++) {
	pmmp2=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m); /*using recursion relation*/
	pmm=pmmp1;
	pmmp1=pmmp2;
      }
      return pmmp2;
    }
  }
}

/*Calculate the theta part of the spherical harmonics (ignoring exp(i*m*phi)).*/
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
