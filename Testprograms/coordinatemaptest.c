/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/poisson.h coordinatemaptest.c */


/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

double boundary(int z, double theta, double phi);
double field(int z, int nr, int nt, int np, double xi, double theta, double phi);
double fieldphysical(int z, double r, double theta, double phi);
double T_n(int n, double x);

int main (void)
{
  FILE *fpgrid;
  FILE *fpfunction;

  int z, i, j, k;
  double xi_i, theta_j, phi_k;
  int nz = 3;
  int nr = 11; /* must be odd? */
  int nt;
  int np = 12; /* must be even */
  int npc;
  scalar2d *b_zjk;
  bound_coeff *bcoeff_zjk;
  scalar2d *f;
  scalar2d *g;
  scalar3d *c_zijk;
  scalar3d *rgrid;
  /*coeff *coeff_zijk;*/
  gsl_vector *alphalist;
  gsl_vector *betalist;
  double r;
  double alpha, beta, bound, bound_map;
  double error;

  nt = np/2 + 1;

  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  b_zjk = scalar2d_alloc(nz-1, nt, np);
  f = scalar2d_alloc(nz, nt, np);
  g = scalar2d_alloc(nz, nt, np);
  bcoeff_zjk = bound_coeff_alloc(nz-1, nt, np);

  c_zijk = scalar3d_alloc(nz, nr, nt, np);
  /*coeff_zijk = coeff_alloc(nz, nr, nt, np);*/

  alphalist = gsl_vector_calloc(nz);
  betalist = gsl_vector_calloc(nz);
 
  rgrid = scalar3d_alloc(nz, nr, nt, np);

  /* Evaluate analytic boundary function on the boundary grid points. */
  for ( z = 0; z < nz-1; z++ ) {
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	scalar2d_set(b_zjk, z, j, k, boundary(z, theta_j, phi_k));
      }
    }
  }

/*   for ( z = 0; z < nz-1; z++ ) { */
/*     for ( j = 0; j < nt; j++ ) { */
/*       theta_j = PI*j/(nt-1); */
/*       for ( k = 0; k < np; k++ ) { */
/* 	phi_k = 2*PI*k/np; */
/* 	printf("(z=%d, j=%d, k=%d, bound=%f, theta=%f, phi=%f)\n", z, j, k, scalar2d_get(b_zjk, z, j, k), theta_j, phi_k); */
/*       } */
/*     } */
/*   } */


  printf("  The boundary b_zjk is:\n");
  print_scalar2d(b_zjk);

  gridtofourier_bound(bcoeff_zjk, b_zjk);
  printf("  The coefficients for the boundary are:\n");
  print_bound_coeff(bcoeff_zjk);

  /* Take boundary points and find functions used for the */
  /* mapping from spherical to surface-matched coordinates. */
  printf("alpha=:\n");
  print_vector(alphalist);
  printf("beta=:\n");
  print_vector(betalist);
  map_physicaltogrid(b_zjk, alphalist, betalist, f, g);
  printf("  f is:\n");
  print_scalar2d(f);
  printf("  g is:\n");
  print_scalar2d(g);
  
  /* kernel bound */
  z = 0;
  for ( j = 0; j < nt; j++ ) {
    theta_j = PI*j/(nt-1);
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      alpha = gsl_vector_get(alphalist, z);
      bound = scalar2d_get(b_zjk, z, j, k);
      bound_map = alpha*(1.0 + scalar2d_get(f, z, j, k) + scalar2d_get(g, z, j, k));
      error = (bound_map - bound) / bound;
      printf("z=%d, j=%d, k=%d, theta=%f, phi=%f, bound_kernel=%f, bound_map_kernel=%f, error=%f\n",
	     z, j, k, theta_j, phi_k, bound, bound_map, error);
    }
  }

  /* inner bound */
  for ( z = 1; z < nz-1; z++ ) {
    for ( j = 0; j < nt; j++ ) {
      theta_j = PI*j/(nt-1);
      for ( k = 0; k < np; k++ ) {
	phi_k = 2*PI*k/np;
	alpha = gsl_vector_get(alphalist, z);
	beta = gsl_vector_get(betalist, z);
	bound = scalar2d_get(b_zjk, z-1, j, k);
	bound_map = alpha*(-1.0 + scalar2d_get(f, z, j, k)) + beta;
	error = (bound_map - bound) / bound;
	printf("z=%d, j=%d, k=%d, theta=%f, phi=%f, bound_in=%f, bound_map_in=%f, error=%f\n",
	       z, j, k, theta_j, phi_k, bound, bound_map, error);
      }
    }
  }
  
  /* outer bound */
  for ( z = 1; z < nz-1; z++ ) {
    for ( j = 0; j < nt; j++ ) {
      theta_j = PI*j/(nt-1);
      for ( k = 0; k < np; k++ ) {
	phi_k = 2*PI*k/np;
	alpha = gsl_vector_get(alphalist, z);
	beta = gsl_vector_get(betalist, z);
	bound = scalar2d_get(b_zjk, z, j, k);
	bound_map = alpha*(1.0 + scalar2d_get(g, z, j, k)) + beta;
	error = (bound_map - bound) / bound;
	printf("z=%d, j=%d, k=%d, theta=%f, phi=%f, bound_out=%f, bound_map_out=%f, error=%f\n",
	       z, j, k, theta_j, phi_k, bound, bound_map, error);
      }
    }
  }
  
  /* external bound */
  z = nz-1;
  for ( j = 0; j < nt; j++ ) {
    theta_j = PI*j/(nt-1);
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      alpha = gsl_vector_get(alphalist, z);
      bound = scalar2d_get(b_zjk, z-1, j, k);
      bound_map = alpha*(-2.0 + scalar2d_get(f, z, j, k));
      error = (bound_map - 1.0/bound) / (1.0/bound);
      printf("z=%d, j=%d, k=%d, theta=%f, phi=%f, 1/bound_kernel=%f, bound_map_kernel=%f, error=%f\n",
	     z, j, k, theta_j, phi_k, 1.0/bound, bound_map, error);
    }
  }
  

/*   /\* Find the radial position of each point. *\/ */
/*   rofxtp(rgrid, alphalist, betalist, f, g); */
/*   printf("  Radial values are:\n"); */
/*   print_scalar3d(rgrid); */
  
/*   /\* print gridpoints in spherical coordinates *\/ */
/*   fpgrid=fopen("gridpoints.txt", "w"); */
/*   for ( z = 0; z < nz-1; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	theta_j = PI*j/(nt-1); */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  phi_k = 2*PI*k/np; */
/* 	  printf("(z=%d, i=%d, j=%d, k=%d, r=%f, theta=%f, phi=%f bound=%f)\n", z, i, j, k, scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k, boundary(z, theta_j, phi_k)); */
/* 	  /\*printf("r=%f\ttheta=%f\tphi=%f\n", scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k);*\/ */
/* 	  /\*fprintf(fpgrid, "%f\t%f\t%f\n", scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k);*\/ */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   /\* invert for external domain *\/ */
/*   z=nz-1; */
/*   for ( i = 0; i < nr-1; i++ ) { */
/*     for ( j = 0; j < nt; j++ ) { */
/*       theta_j = PI*j/(nt-1); */
/*       for ( k = 0; k < np; k++ ) { */
/* 	phi_k = 2*PI*k/np; */
/* 	/\*printf("r=%f\ttheta=%f\tphi=%f\n", scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k); */
/* 	/\*fprintf(fpgrid, "%f\t%f\t%f\n", 1/scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k);*\/ */
/*       } */
/*     } */
/*   } */
  
  
/*   /\* Assign the function data. *\/ */
/*   for ( z = 0; z < nz-1; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  if(z==0) */
/* 	    xi_i = sin(PI*i/(2*(nr-1))); */
/* 	  else */
/* 	    xi_i = -cos(PI*i/(nr-1)); */
/* 	  theta_j = PI*j/(nt-1); */
/* 	  phi_k = 2*PI*k/np; */
/* 	  r = scalar3d_get(rgrid, z, i, j, k); */
/* 	  scalar3d_set(c_zijk, z, i, j, k, fieldphysical(z, r, theta_j, phi_k)); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   z=nz-1; */
/*   for ( i = 0; i < nr-1; i++ ) { */
/*     for ( j = 0; j < nt; j++ ) { */
/*       for ( k = 0; k < np; k++ ) { */
/* 	xi_i = -cos(PI*i/(nr-1)); */
/* 	theta_j = PI*j/(nt-1); */
/* 	phi_k = 2*PI*k/np; */
/* 	r = scalar3d_get(rgrid, z, i, j, k); */
/* 	scalar3d_set(c_zijk, z, i, j, k, fieldphysical(z, 1/r, theta_j, phi_k)); */
/*       } */
/*     } */
/*   } */
/*   z=nz-1; */
/*   i=nr-1; */
/*   for ( j = 0; j < nt; j++ ) { */
/*     for ( k = 0; k < np; k++ ) { */
/*       scalar3d_set(c_zijk, z, i, j, k, 0); */
/*     } */
/*   } */
  
/*   /\*print_scalar3d(c_zijk);*\/ */
  
/*   /\* print to file the function and radial position for all gridpoints *\/ */
/*   fpgrid=fopen("fofr.txt", "w"); */
/*   for ( z = 0; z < nz-1; z++ ) { */
/*     for ( i = 0; i < nr; i++ ) { */
/*       for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  fprintf(fpgrid, "%f\t%f\n", scalar3d_get(rgrid, z, i, j, k), scalar3d_get(c_zijk, z, i, j, k)); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   z=nz-1; */
/*   for ( i = 0; i < nr-1; i++ ) { */
/*     for ( j = 0; j < nt; j++ ) { */
/*       for ( k = 0; k < np; k++ ) { */
/* 	fprintf(fpgrid, "%f\t%f\n", 1/scalar3d_get(rgrid, z, i, j, k), scalar3d_get(c_zijk, z, i, j, k)); */
/*       } */
/*     } */
/*   } */

  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.1*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else */
/*     return 5.0; */
/* } */

double boundary(int z, double theta, double phi)
{
  if(z==0)
    return 2.0*(1.0 + /*0.3*sin(theta)*(cos(phi)+sin(phi)) +*/ 0.2*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi)));
  else
    return 5.0;
}

/* double boundary(int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 2.0; */
/*   else */
/*     return 5.0; */
/* } */

/* double boundary(int nt, int np, double theta, double phi) */
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
  return cos(r);
}


/******************************************************/
/* Evaluate a function at the point (xi, theta, phi). */
/******************************************************/
double field(int z, int nr, int nt, int np, double xi, double theta, double phi)
{
  int i, j, k;
  double basis; /* basis function */
  double sum;
  int npc = np/2 + 1;
  
  sum = 0.0;
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < npc; k++ ) {
      
      if( z==0 && j%2 ) { /* odd j, zone 0 */
	for ( i = 0; i < nr-1; i++ ) { /* avoid aliasing (nr not allowed) */
	  
	  basis = cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi);
	  
	  if(!(k%2)) { /* even k */
	    basis *= cos(j*theta);
	  } else { /* odd k */
	    basis *= sin(j*theta);
	  }
	  
	  basis *= T_n(2*i+1, xi);
	  
	  sum += 1.1111111*basis;
	}
      } else if( z==0 && !(j%2) ) { /* odd j, zone 0 */
	for ( i = 0; i < nr; i++ ) {
	  
	  basis = cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi);
	  
	  if(!(k%2)) { /* even k */
	    basis *= cos(j*theta);
	  } else { /* odd k */
	    basis *= sin(j*theta);
	  }
	  
	  basis *= T_n(2*i, xi);
	  
	  sum += 1.1111111*basis;
	}
      } else { /* other zones */
	for ( i = 0; i < nr; i++ ) {
	  
	  basis = cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi);
	  
	  if(!(k%2)) { /* even k */
	    basis *= cos(j*theta);
	  } else { /* odd k */
	    basis *= sin(j*theta);
	  }
	  
	  basis *= T_n(i, xi);
	  
	  sum += 1.1111111*basis;
	}
      }
    }
  }
  return sum;
  
/*   sum = 0.0; */
/*   for ( i = 0; i < nr-1; i++ ) { */
  
/* 	basis = cos(0*phi) + sin(0*phi); */
  
/* 	basis *= cos(1*theta); */
  
  
  
/* 	basis *= T_n(2*i+1, xi); */
  
/* 	sum += 1.1111111*basis; */
  
/*   } */
/*   return sum; */
  
  /*return 1.1111111*T_n(2, xi)*cos(2*theta)*(cos(2*phi) + sin(2*phi));*/
}


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
