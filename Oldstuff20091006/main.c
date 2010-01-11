/********************************************************************/
/* TO DO: 1) chebfit_odd, chebeval_odd functions are not working    */
/*        2) use LAPACK to solve pentadiagonal matrices             */
/********************************************************************/

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "poisson.h"

/* function prototypes */
double boundary(int zone, double theta, double phi);
double source(int zone, double r, double theta, double phi);

/* main program */
int main (void)
{
  grid3d *points = grid3d_alloc(2, 5, 1, 1);
  scalar2d *b_zjk = scalar2d_alloc(1, 1, 1);
  bound_coeff *b_zlm = bound_coeff_alloc(1, 1);
  scalar3d *s_zijk = scalar3d_alloc(2, 5, 1, 1);
  coeff *c_znlm = coeff_alloc(2, 5, 1);
  coeff *f_znlm = coeff_alloc(2, 5, 1);

  makegrid3d(points);
  /*print_grid3d(points);*/
  
  boundtogrid(b_zjk, points, boundary);
  /*print_scalar2d(b_zjk);*/
  
  decompose_bound(b_zlm, b_zjk);
  printf("boundary coefficients\n");
  print_bound_coeff(b_zlm);
  
  ftogrid(s_zijk, b_zjk, points, source);
  /*print_scalar3d(s_zijk);*/
  
  decompose(c_znlm, s_zijk);
  printf("source coefficients\n");
  print_coeff_2(c_znlm);
  
  solve_poisson(f_znlm, c_znlm, b_zlm);
  printf("solution coefficients\n");
  print_coeff_2(f_znlm);  
  
  return 0;
}


double boundary(int zone, double theta, double phi)
{
  return 1;

  /*if(zone==0)
	return 1; 
  else if(zone==1)
	return 2;
  else
  return 4;*/
  
  /*  if(zone==0)
	return 1*sqrt(1/(4*PI)); 
  else if(zone==1)
	return 17*sqrt(1/(4*PI)) 
	  + 5*sqrt(3/(4*PI))*cos(theta) 
	  + 7*(-1)*sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi)
	  + 8*0.25*sqrt(105/(2*PI))*sin(theta)*sin(theta)*cos(theta)*sin(2*phi);
  else
	return 20*sqrt(1/(4*PI)); */
}


double source(int zone, double r, double theta, double phi)
{
  if(zone==0)
	return 4*PI*1;
  else
	return 0;
  
  /*if(zone==0)
	return 1*sqrt(1/(4*PI))
	+ 5*(r)*sqrt(3/(4*PI))*cos(theta)
	+ 6*(2*r*r-1)*(-1)*sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi);
	else
	return 0;*/
  
  /*if(zone==0)
	return 1; 
	else if(zone==1)
	return 1;
	else if(zone==2)
	return 1;
  else
  return 1;*/
  
/*   if(zone==0) */
/* 	return 1*sqrt(1/(4*PI));  */
/*   else if(zone==1) */
/* 	return (2*r*r-1)*(1*sqrt(1/(4*PI)))  */
/* 	  + 2*sqrt(3/(4*PI))*cos(theta)  */
/* 	  + 3*(-1)*sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi); */
/*   else if(zone==2) */
/* 	return 4*sqrt(1/(4*PI))  */
/* 	  + 5*(2*r*r-1)*sqrt(3/(4*PI))*cos(theta)  */
/* 	  + 6*(-1)*sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi); */
/*   else */
/* 	return 7*sqrt(1/(4*PI))  */
/* 	  + 8*sqrt(3/(4*PI))*cos(theta)  */
/* 	  + 9*(4*r*r*r-3*r)*(-1)*sqrt(15/(8*PI))*sin(theta)*cos(theta)*sin(phi); */
}
