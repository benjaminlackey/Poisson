/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/poisson.h dfdthetatest.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

double boundary(int z, int nt, int np, double theta, double phi);
double field(int z, int nr, int nt, int np, double xi, double theta, double phi);
double fieldphysical(int z, double r, double theta, double phi);
double gradfield_theta(int z, double r, double theta, double phi);
double T_n(int n, double x);

int main (void)
{

 
  /*clock_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12; 
    float ratio;*/ 

  FILE *fpgrid;
  /*FILE *fpfunction;*/

  int z, i, j, k;
  int nz = 3;
  int nr = 5; /* must be odd? */
  int nt = 7;
  int np = 12; /* must be even */
  int npc;
  double xi_i, r_i, theta_j, phi_k;
  double r;
  scalar2d *surface_grid;
  scalar2d *f_grid;
  scalar2d *g_grid;
  gsl_vector *alphalist;
  gsl_vector *betalist;
  scalar3d *r_grid;
  scalar3d *field_grid;
  coeff *field_coeff;
  coeff *dfdxi_coeff;
  scalar3d *dfdxi_grid;
  scalar3d *j1_grid;
  scalar3d *j2_grid;
  scalar3d *temp_grid;
  scalar3d *fielderase_grid;
  coeff *dfdt_coeff;
  scalar3d *onebyrdfdt_grid;
  scalar3d *r_scalar3d;

  double grad, grad_analytic, error;

  npc = ( np / 2 ) + 1;

  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  surface_grid = scalar2d_alloc(nz-1, nt, np);
  f_grid = scalar2d_alloc(nz, nt, np);
  g_grid = scalar2d_alloc(nz, nt, np);
  alphalist = gsl_vector_calloc(nz);
  betalist = gsl_vector_calloc(nz);
  
  /* value of R (or U) at each point (xi, theta, phi) on grid */
  r_grid = scalar3d_alloc(nz, nr, nt, np);
  
  /* the field */
  field_grid = scalar3d_alloc(nz, nr, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);
  
  /* derivatives of field */
  dfdxi_grid = scalar3d_alloc(nz, nr, nt, np);
  dfdxi_coeff = coeff_alloc(nz, nr, nt, np); 
  dfdt_coeff = coeff_alloc(nz, nr, nt, np);

  j1_grid = scalar3d_alloc(nz, nr, nt, np);
  j2_grid = scalar3d_alloc(nz, nr, nt, np);
  temp_grid = scalar3d_alloc(nz, nr, nt, np); 
  fielderase_grid = scalar3d_alloc(nz, nr, nt, np);
  onebyrdfdt_grid = scalar3d_alloc(nz, nr, nt, np);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);

  /* Evaluate analytic boundary function on the boundary grid points. */
  for ( z = 0; z < nz-1; z++ ) {
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	scalar2d_set(surface_grid, z, j, k, boundary(nt, np, z, theta_j, phi_k));
      }
    }
  }
  /*print_scalar2d(surface_grid);*/

  /* determine boundary functions given the boundary surfaces */
  map_physicaltogrid_kernel(surface_grid, alphalist, f_grid, g_grid);
  map_physicaltogrid_shell(surface_grid, 1, alphalist, betalist, f_grid, g_grid);
  map_physicaltogrid_ext(surface_grid, alphalist, f_grid);
/*   print_vector(alphalist); */
/*   print_vector(betalist); */
/*   print_scalar2d(f_grid); */
/*   print_scalar2d(g_grid); */
  /* Find the radial position of each point. */
  rofxtp(r_grid, alphalist, betalist, f_grid, g_grid);

  /* Assign the field data to field_grid. */
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if(z==0)
	    xi_i = sin(PI*i/(2*(nr-1)));
	  else
	    xi_i = -cos(PI*i/(nr-1));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  r = scalar3d_get(r_grid, z, i, j, k);
	  scalar3d_set(field_grid, z, i, j, k, fieldphysical(z, r, theta_j, phi_k));
	}
      }
    }
  }
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	xi_i = -cos(PI*i/(nr-1));
	theta_j = PI*j/(nt-1);
	phi_k = 2*PI*k/np;
	r = scalar3d_get(r_grid, z, i, j, k);
	scalar3d_set(field_grid, z, i, j, k, fieldphysical(z, 1/r, theta_j, phi_k));
      }
    }
  }
  z=nz-1;
  i=nr-1;
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      scalar3d_set(field_grid, z, i, j, k, 0);
    }
  }
  /*print_scalar3d(field_grid);*/
  jacobian1(j1_grid, alphalist, f_grid, g_grid);
  /*print_scalar3d(j1_grid);*/
  jacobian2(j2_grid, alphalist, betalist, f_grid, g_grid);
  /*print_scalar3d(j1_grid);*/

  /* transforms erase grid so set it to another erasable grid */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(fielderase_grid, z, i, j, k, scalar3d_get(field_grid, z, i, j, k));
	}
      }
    }
  }  
  /* evaluate df/dxi */
  gridtofourier(field_coeff, fielderase_grid, 0, 0);
  dfdxi(dfdxi_coeff, field_coeff);
  fouriertogrid(dfdxi_grid, dfdxi_coeff, 1, 0);
  

  /* transforms erase grid so set it to another erasable grid */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(fielderase_grid, z, i, j, k, scalar3d_get(field_grid, z, i, j, k));
	}
      }
    }
  }  
  /* evaluate 1/R * df/dtheta */
  gridtofourier(field_coeff, fielderase_grid, 0, 0);
  dfdthetaprime(dfdt_coeff, field_coeff);
  dividebyr(onebyrdfdt_grid, dfdt_coeff, 0, 1, alphalist, betalist, f_grid, g_grid);
  
/*   printf("field_grid:\n"); */
/*   print_scalar3d(field_grid); */
/*   printf("j2_grid:\n");  */
/*   print_scalar3d(j2_grid);  */
/*   printf("onebyrdfdt_grid:\n"); */
/*   print_scalar3d(onebyrdfdt_grid);  */
  
  /* set 1/r * df/dtheta = (1/R)*dfdt' - (J_2/J_1)*df/dxi */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(temp_grid, z, i, j, k, scalar3d_get(onebyrdfdt_grid, z, i, j, k) - scalar3d_get(j2_grid, z, i, j, k)*scalar3d_get(dfdxi_grid, z, i, j, k)/scalar3d_get(j1_grid, z, i, j, k));
	}
      }
    }
  }

  /*print_scalar3d(temp_grid);*/
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alphalist, betalist, f_grid, g_grid);
  print_scalar3d(r_scalar3d);
/*   for ( z = 0; z < nz-1; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  r_i = ((z==nz-1) ? 1.0/scalar3d_get(r_scalar3d, z, i, j, k) : scalar3d_get(r_scalar3d, z, i, j, k)); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/* 	  grad = scalar3d_get(temp_grid, z, i, j, k); */
/* 	  grad_analytic = gradfield_theta(z, r_i, theta_j, phi_k); */
/* 	  error = (grad - grad_analytic)/grad_analytic; */
/* 	  printf("z=%d, i=%d, j=%d, k=%d, r_i=%.18e, t_j=%.18e, p_k=%.18e, %.18e, %.18e, %.18e\n", z, i, j, k, r_i, theta_j, phi_k, grad, grad_analytic, error); */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   /\* print to file the function and radial position for all gridpoints *\/ */
/*   fpgrid=fopen("onebyrdfdtheta.txt", "w"); */
/*   for ( z = 0; z < nz-1; z++ ) { */
/*     for ( j = 0; j < nt; j++ ) { */
/*       theta_j = PI*j/(nt-1); */
/*       for ( i = 1; i < nr; i++ ) { /\* dfdtheta is not well defined for r = 0? *\/  */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  fprintf(fpgrid, "%.18e\t%.18e\n", theta_j, scalar3d_get(temp_grid, z, i, j, k)); */
/* 	} */
/*       } */
/*     } */
/*   } */
  
  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0; */
/*   else */
/*     return 5.0; */
/* } */

double boundary(int nt, int np, int z, double theta, double phi)
{
  if(z==0)
    return 1.0*(1.0 + 0.3*sin(theta)*(cos(phi)+sin(phi)) + 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi)));
  else
    return 5.0*(1.0 - 0.2*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi)));
}

/* double boundary(int z, int nt, int np, double theta, double phi) */
/* { */
/*   int j, k; */
/*   double sum; */
/*   int npc = np/2 + 1; */

/*   sum = 1; */
/*   for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < npc; k+=2 ) { */
/* 	  sum += /\*((j+1)*10+(k+1))*\/1.1111111*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/* 	} */
/* 	for ( k = 1; k < npc; k+=2 ) { */
/* 	  sum += /\*((j+1)*10+(k+1))*\/1.1111111*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi)); */
/* 	} */
/*   } */
  
/*   /\*sum = 0; */
/*   for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np/2+1; k+=2 ) { */
/* 	  sum += (j+1)*cos(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	} */
/* 	for ( k = 1; k < np/2+1; k+=2 ) { */
/* 	  sum += (j+1)*sin(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	} */
/* 	}*\/ */
  
/*   /\*sum = 0; */
/* 	for ( j = 0; j < nt; j++ ) { */
/* 	sum += (j+1)*cos(j*theta); */
/* 	}*\/ */
  
/*   return sum; */
/* } */


/*********************************************************/
/* Some function of position in spherical coordinates.   */
/*********************************************************/
double fieldphysical(int z, double r, double theta, double phi)
{
  return r*cos(theta);
}

double gradfield_theta(int z, double r, double theta, double phi)
{
    return -sin(theta);
}

/******************************************************/
/* Evaluate a function at the point (xi, theta, phi). */
/******************************************************/
double field(int z, int nr, int nt, int np, double xi, double theta, double phi)
{
  return cos(xi);
}

/* double field(int z, int nr, int nt, int np, double xi, double theta, double phi) */
/* { */
/*   int i; */
/*   double sum; */
  
/*   if(z==0) { */
/*     sum = 0.0; */
/*     for(i=0; i<nr; i++){ */
/*       sum += T_n(2*i, xi); */
/*     } */
/*     for(i=0; i<nr-1; i++){ */
/*       sum += T_n(2*i+1, xi)*(cos(theta) + sin(theta)*(cos(phi)+sin(phi))); */
/*     } */
/*     return sum; */
/*   } else { */
/*     sum = 0.0; */
/*     for(i=0; i<nr; i++){ */
/*       sum += T_n(i, xi); */
/*     } */
/*     return sum; */
/*   } */
/* } */


/************************************************/
/* Evaluate nth Chebyshev polynomial at point x */
/************************************************/
double T_n(int n, double x)
{
  double Tn;          /* nth Chebyshev polynomial T_n(x) */
  double Tnminus1;    /* (n-1)th Chebyshev polynomial T_{n-1}(x) */ 
  double Tnminus2;    /* (n-2)th Chebyshev polynomial T_{n-2}(x) */ 
  int j; 
    
  /* Evaluate nth Chebyshev polynomial at point x.     */
  /* Use the recursion relation:                       */
  /* T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x)              */

  if(n==0) {
	return 1;
  } else if(n==1) {
	return x;
  } else {
	Tnminus2 = 1;
	Tnminus1 = x;
	for(j=2; j<=n; j++) {
	  Tn = 2.0*x*Tnminus1 - Tnminus2;
	  Tnminus2 = Tnminus1;
	  Tnminus1 = Tn;
	}
	return Tn;
  } 
}
