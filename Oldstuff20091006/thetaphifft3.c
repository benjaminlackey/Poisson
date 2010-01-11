/* To compile type: gcc -I/opt/local/include -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas thetaphifft3.c print.c coefficients.c poisson.h */


/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* fftw header */
#include <fftw3.h>

#include "poisson.h"

#define PI 3.141592653589793

#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
/* #define iseven(i) ((i)%2 ? 0 : 1) /\* 0 if i is odd. 1 if i is even *\/ */
/* #define isodd(i) ((i)%2 ? 1 : 0) /\* 1 if i is odd. 0 if i is even *\/ */

double boundary(int nt, int np, int z, double theta, double phi);
void gridtofourier_bound(bound_coeff *bcoeff_zjk, scalar2d *b_zjk);
void fouriertogrid_bound(scalar2d *b_zjk, bound_coeff *bcoeff);
void gridtofourier(coeff *coeff, scalar3d *c_zijk);
void fouriertogrid(scalar3d *c_zijk, coeff *coeff);
void map_physicaltogrid_kernel(scalar2d *b_zjk, gsl_vector *alphalist, scalar2d *fodd_grid, scalar2d *geven_grid);
void map_physicaltogrid_shell(scalar2d *b_zjk, int z, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *fin_grid, scalar2d *gout_grid);
void map_physicaltogrid_ext(scalar2d *b_zjk, gsl_vector *alphalist, scalar2d *f_grid);
void rofxtp(scalar3d *rgrid, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g);
scalar3d jacobian1(scalar2d *f, scalar2d *g);
scalar3d jacobian2(scalar2d *f, scalar2d *g);
scalar3d jacobian3(scalar2d *f, scalar2d *g);
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
  int nt = 7;
  int np = 6; /* must be even */
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
	scalar2d_set(b_zjk, z, j, k, boundary(nt, np, z, theta_j, phi_k));
      }
    }
  }

  printf("  The boundary b_zjk is:\n");
  print_scalar2d(b_zjk);

  gridtofourier_bound(bcoeff_zjk, b_zjk);
  printf("  The coefficients for the boundary are:\n");
  print_bound_coeff(bcoeff_zjk);

  /* Take boundary points and find functions used for the */
  /* mapping from spherical to surface-matched coordinates. */
/*   printf("  For kernel, fodd_grid starts off as:\n"); */
/*   print_scalar2d(f); */
/*   printf("  For kernel, geven_grid starts off as:\n"); */
/*   print_scalar2d(g); */
/*   print_vector(alphalist); */
/*   print_vector(betalist); */
  map_physicaltogrid_kernel(b_zjk, alphalist, f, g);
  printf("  f is:\n");
  print_scalar2d(f);
  printf("  g is:\n");
  print_scalar2d(g);

  map_physicaltogrid_shell(b_zjk, 1, alphalist, betalist, f, g);
  printf("  f is:\n");
  print_scalar2d(f);
  printf("  g is:\n");
  print_scalar2d(g);
  
  map_physicaltogrid_ext(b_zjk, alphalist, f);
  printf("  f is:\n");
  print_scalar2d(f);
  printf("  g is:\n");
  print_scalar2d(g);
  
  /* Find the radial position of each point. */
  rofxtp(rgrid, alphalist, betalist, f, g);
  printf("  Radial values are:\n");
  print_scalar3d(rgrid);
  
  /* print gridpoints in spherical coordinates */
  fpgrid=fopen("gridpoints.txt", "w");
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	for ( k = 0; k < np; k++ ) {
	  phi_k = 2*PI*k/np;
	  /*printf("r=%f\ttheta=%f\tphi=%f\n", scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k);*/
	  fprintf(fpgrid, "%f\t%f\t%f\n", scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k);
	}
      }
    }
  }
  /* invert for external domain */
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      theta_j = PI*j/(nt-1);
      for ( k = 0; k < np; k++ ) {
	phi_k = 2*PI*k/np;
	/*printf("r=%f\ttheta=%f\tphi=%f\n", scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k);*/
	fprintf(fpgrid, "%f\t%f\t%f\n", 1/scalar3d_get(rgrid, z, i, j, k), theta_j, phi_k);
      }
    }
  }
  
  
  /* Assign the function data. */
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
	  r = scalar3d_get(rgrid, z, i, j, k);
	  scalar3d_set(c_zijk, z, i, j, k, fieldphysical(z, r, theta_j, phi_k));
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
	r = scalar3d_get(rgrid, z, i, j, k);
	scalar3d_set(c_zijk, z, i, j, k, fieldphysical(z, 1/r, theta_j, phi_k));
      }
    }
  }
  z=nz-1;
  i=nr-1;
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      scalar3d_set(c_zijk, z, i, j, k, 0);
    }
  }
  
  print_scalar3d(c_zijk);
  
  /* print to file the function and radial position for all gridpoints */
  fpgrid=fopen("fofr.txt", "w");
  for ( z = 0; z < nz-1; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  fprintf(fpgrid, "%f\t%f\n", scalar3d_get(rgrid, z, i, j, k), scalar3d_get(c_zijk, z, i, j, k));
	}
      }
    }
  }
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	fprintf(fpgrid, "%f\t%f\n", 1/scalar3d_get(rgrid, z, i, j, k), scalar3d_get(c_zijk, z, i, j, k));
      }
    }
  }



/*   npc = ( np / 2 ) + 1; */
  
/*   /\* Assign the boundary data. *\/ */
/*   for ( z = 0; z < nz-1; z++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  phi_k = 2*PI*k/np; */
/* 	  for ( j = 0; j < nt; j++ ) { */
/* 		theta_j = PI*j/(nt-1); */
/* 		scalar2d_set(b_zjk, z, j, k, boundary(nt, np, theta_j, phi_k)); */
/* 	  } */
/* 	} */
/*   } */
  
/* /\*   /\\* Assign the boundary data version 2. *\\/ *\/ */
/* /\*   for ( z = 0; z < nz-1; z++ ) { *\/ */
/* /\* 	for ( k = 0; k < np; k++ ) { *\/ */
/* /\* 	  scalar2d_set(b_zjk, z, 0, k, (z+1)*100); *\/ */
/* /\* 	  for ( j = 1; j < nt-1; j++ ) { *\/ */
/* /\* 		scalar2d_set(b_zjk, z, j, k, (z+1)*100 + (j+1)*10 + (k+1)); *\/ */
/* /\* 	  } *\/ */
/* /\* 	  scalar2d_set(b_zjk, z, nt-1, k, (z+1)*100); *\/ */
/* /\* 	} *\/ */
/* /\*   } *\/ */
  
/*   /\* Assign the field data. *\/ */
/*   for ( z = 0; z < nz; z++ ) { */
/* 	for ( i = 0; i < nr; i++ ) { */
/* 	  for ( j = 0; j < nt; j++ ) { */
/* 		for ( k = 0; k < np; k++ ) { */
/* 		  if(z==0) { */
/* 			xi_i = sin(PI*i/(2*(nr-1))); */
/* 		  } else { */
/* 			xi_i = -cos(PI*i/(nr-1)); */
/* 		  } */
/* 		  theta_j = PI*j/(nt-1);  */
/* 		  phi_k = 2*PI*k/np; */
/* 		  scalar3d_set(c_zijk, z, i, j, k, field(z, nr, nt, np, xi_i, theta_j, phi_k)); */
/* 		} */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* print start data*\/ */
/*   /\*printf("  Boundary data.\n");*\/ */
  
/*   /\*print_scalar2d(b_zjk);*\/ */
  
/*   gridtofourier_bound(bcoeff_zjk, b_zjk); */
  
/*   /\*print_bound_coeff(bcoeff_zjk);*\/ */
  
/*   fouriertogrid_bound(b_zjk, bcoeff_zjk); */

/*   /\*print_scalar2d(b_zjk);*\/ */

/*   /\*printf("  Field data.\n");*\/ */

/*   /\*print_scalar3d(c_zijk);*\/ */

/*   gridtofourier(coeff_zijk, c_zijk); */

/*   print_coeff(coeff_zijk); */

/*   fouriertogrid(c_zijk, coeff_zijk); */
  
/*   print_scalar3d(c_zijk); */
  
/*   gridtofourier(coeff_zijk, c_zijk);  */
  
/*   print_coeff(coeff_zijk); */
  
  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
double boundary(int nt, int np, int z, double theta, double phi)
{
  if(z==0)
    return 1.0*(1.0 + 0.1*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi)));
  else
    return 5.0;
}
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


/*************************************************************************/
/* Find the scale factor alpha and functions fodd, geven                 */
/* used in the mapping from physical points to collocation grid points.  */
/* This is for the kernel.                                               */
/*************************************************************************/
void map_physicaltogrid_kernel(scalar2d *b_zjk, gsl_vector *alphalist, scalar2d *fodd_grid, scalar2d *geven_grid)
{
  int z;
  int j;
  int k;
  int nz;
  int nt;
  int np;
  int npc;
  double specialpoint;
  double mu;
  double new;
  double alpha;
  scalar2d *f_grid; /* rescaled boundary */
  bound_coeff *f_fourier; /* rescaled boundary in fourier space */
  bound_coeff *fodd_fourier; /* odd phi components (e^{i*k*phi}) */
  bound_coeff *geven_fourier; /* even phi components (e^{i*k*phi}) */
  
  /* nz for b_zjk only has nz-1 values */
  nz = b_zjk->nz + 1;
  nt = b_zjk->nt;
  np = b_zjk->np;
  npc = np/2 + 1;
  
  /*if(np%2 == 1)
	printf("np must be even in gridtofourier_bound\n");*/
  
  /* there are nz-1 boundaries */
  f_grid = scalar2d_alloc(nz-1, nt, np);
  f_fourier = bound_coeff_alloc(nz-1, nt, np);
  /* there are nz zones each of which have functions f, g */
  fodd_fourier = bound_coeff_alloc(nz, nt, np);
  geven_fourier = bound_coeff_alloc(nz, nt, np);
/*   printf("fodd_fourier starts off as:\n"); */
/*   print_bound_coeff(fodd_fourier); */
/*   printf("geven_fourier starts off as:\n"); */
/*   print_bound_coeff(geven_fourier); */
  
  /* Initialize these structures to 0 */
  /* you should probably just make an initialization function instead */
  for(z=0; z<nz-1; z++) {
    for(j=0; j<nt; j++) {
      for(k=0; k<npc; k++) {
	/*printf("z=%d, j=%d, k=%d\n", z, j, k);*/
	scalar2d_set(f_grid, z, j, k, 0.0);
	bound_coeff_set(f_fourier, z, j, k, REAL, 0.0);
	bound_coeff_set(f_fourier, z, j, k, IMAG, 0.0);
      }
    }
  }
  for(z=0; z<nz; z++) {
    for(j=0; j<nt; j++) {
      for(k=0; k<npc; k++) {
	/*printf("z=%d, j=%d, k=%d\n", z, j, k);*/
	bound_coeff_set(fodd_fourier, z, j, k, REAL, 0.0);
	bound_coeff_set(fodd_fourier, z, j, k, IMAG, 0.0);
	bound_coeff_set(geven_fourier, z, j, k, REAL, 0.0);
	bound_coeff_set(geven_fourier, z, j, k, IMAG, 0.0);
      }
    }
  }
  
  /* evaluate some special point and find f = S - S_point for each point on kernel boundary */
  specialpoint = scalar2d_get(b_zjk, 0, 0, 0);
  printf("specialpoint=%f\n", specialpoint);
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
	  scalar2d_set( f_grid, 0, j, k, 
			scalar2d_get(b_zjk, 0, j, k) - specialpoint );
    }
  }
  printf("f_grid = b_zjk - specialpoint\n");
  print_scalar2d(f_grid);

  /* decompose f into double fourier series */
  gridtofourier_bound(f_fourier, f_grid);
  print_bound_coeff(f_fourier);

  /* set fodd~ to odd phi harmonics */
  for(k=1; k<npc; k += 2) {
    for(j=0; j<nt; j++) {
      bound_coeff_set(fodd_fourier, 0, j, k, REAL, bound_coeff_get(f_fourier, 0, j, k, REAL));
      bound_coeff_set(fodd_fourier, 0, j, k, IMAG, bound_coeff_get(f_fourier, 0, j, k, IMAG));
    }
  }
  printf("fodd_fourier:\n");
  print_bound_coeff(fodd_fourier);

  /* set geven~ to even phi harmonics */
  for(k=0; k<npc; k += 2) {
    for(j=0; j<nt; j++) {
      bound_coeff_set(geven_fourier, 0, j, k, REAL, bound_coeff_get(f_fourier, 0, j, k, REAL));
      bound_coeff_set(geven_fourier, 0, j, k, IMAG, bound_coeff_get(f_fourier, 0, j, k, IMAG));
    }
  }
  printf("geven_fourier:\n");
  print_bound_coeff(geven_fourier);  

  
  /* transform fodd~ and geven~ back to spatial domain */
  fouriertogrid_bound(fodd_grid, fodd_fourier);
  fouriertogrid_bound(geven_grid, geven_fourier);
  printf("fodd_grid:\n");
  print_scalar2d(fodd_grid);
  printf("fodd_grid:\n");
  print_scalar2d(fodd_grid);
  
  /* calculate mu (the negative of the minimum value of geven~) */
  mu = scalar2d_get(geven_grid, 0, 0, 0);
  for(k=0; k<np; k++) {
    for(j=0; j<nt; j++) {
      new = scalar2d_get(geven_grid, 0, j, k);
      mu = (new < mu ? new : mu);
    }
  }
  mu = -mu;
  printf("mu=%f\n", mu);

  /* calculate alpha */
  alpha = specialpoint - mu;
  gsl_vector_set(alphalist, 0, alpha);
  printf("alpha=%f\n", alpha);

  /* rescale fodd~ and geven~ to calculate fodd and geven */ 
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
      /* fodd: */
      scalar2d_set(fodd_grid, 0, j, k, scalar2d_get(fodd_grid, 0, j, k)/alpha);
      /* geven: */
      scalar2d_set(geven_grid, 0, j, k, scalar2d_get(geven_grid, 0, j, k)/alpha + mu);
    }
  }
  
  /* you still need to free the following things: */
  /*scalar2d_free(f_grid);
  bound_coeff_free(f_fourier);
  bound_coeff_free(fodd_fourier);
  bound_coeff_free(geven_fourier); */
  /*printf("oh hi\n");*/
}


/*************************************************************************/
/* Find the scale factors alpha, beta and functions fin, gout            */
/* used in the mapping from physical points to collocation grid points.  */
/* This is for the shells.                                               */
/*************************************************************************/
void map_physicaltogrid_shell(scalar2d *b_zjk, int z, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *fin_grid, scalar2d *gout_grid)
{
  int j;
  int k;
  int nt;
  int np;
  double specialpointin;
  double specialpointout;
  double lambda;
  double mu;
  double new;
  double alpha;
  double beta;
  
  nt = b_zjk->nt;
  np = b_zjk->np;
  
  /*if(np%2 == 1)
    printf("np must be even in gridtofourier_bound\n");*/
  
  /* evaluate some special point on inner and outer surfaces and find */ 
  specialpointin = scalar2d_get(b_zjk, z-1, 0, 0);
  specialpointin = scalar2d_get(b_zjk, z, 0, 0);
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
      /* f~ = Sin - Sin_point */
      /* fill in zone z of fin_grid even though this corresponds to boundary z-1 */
      scalar2d_set( fin_grid, z, j, k, 
		    scalar2d_get(b_zjk, z-1, j, k) - specialpointin );
      /* g~ = Sout - Sout_point */	
      scalar2d_set( gout_grid, z, j, k, 
		    scalar2d_get(b_zjk, z, j, k) - specialpointout );
    }
  }
  
  /* calculate lambda (the negative of the maximum value of f~) */
  lambda = scalar2d_get(fin_grid, z, 0, 0);
  for(k=0; k<np; k++) {
    for(j=0; j<nt; j++) {
      new = scalar2d_get(fin_grid, z, j, k);
      lambda = (new > lambda ? new : lambda);
    }
  }
  lambda = -lambda;
  
  /* calculate mu (the negative of the minimum value of g~) */
  mu = scalar2d_get(gout_grid, z, 0, 0);
  for(k=0; k<np; k++) {
    for(j=0; j<nt; j++) {
      new = scalar2d_get(gout_grid, z, j, k);
      mu = (new < mu ? new : mu);
    }
  }
  mu = -mu;
  
  /* calculate alpha and beta */
  alpha = 0.5*(specialpointout - specialpointin + lambda - mu);
  gsl_vector_set(alphalist, z, alpha);
  beta = 0.5*(specialpointout + specialpointin - lambda - mu);
  gsl_vector_set(betalist, z, beta);

  /* rescale fin~ and gout~ to calculate fin and gout */ 
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
      /* fin: */
      scalar2d_set(fin_grid, z, j, k, (scalar2d_get(fin_grid, z, j, k) + lambda)/alpha);
      /* gout: */
      scalar2d_set(gout_grid, z, j, k, (scalar2d_get(gout_grid, z, j, k) + mu)/alpha);
    }
  }
  
}


/*************************************************************************/
/* Find the scale factor alpha and function f                            */
/* used in the mapping from physical points to collocation grid points.  */
/* This is for the external domain.                                      */
/*************************************************************************/
void map_physicaltogrid_ext(scalar2d *b_zjk, gsl_vector *alphalist, scalar2d *f_grid)
{
  int j;
  int k;
  int nz; /* nz boundaries implies nz+1 zones */
  int nt;
  int np;
  double specialpointin;
  double lambda;
  double new;
  double alpha;
  
  /* number of zones is 1 more than number of boundaries */
  nz = b_zjk->nz + 1;
  nt = b_zjk->nt;
  np = b_zjk->np;
  
  /*if(np%2 == 1)
	printf("np must be even in gridtofourier_bound\n");*/
  
  /*printf("hi again\n");*/
  /* evaluate some special point on inner surface */ 
  specialpointin = scalar2d_get(b_zjk, nz-2, 0, 0);
  printf("specialpointin=%f\n", specialpointin);
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
      /*printf("j=%d, k=%d\n", j, k);*/
      /* f~ = Sin^-1 - Sin_point^-1 */
      /* fill in zone nz-1 of f_grid even though this is given by boundary nz-2 */
      scalar2d_set( f_grid, nz-1, j, k, 
		    1.0/scalar2d_get(b_zjk, nz-2, j, k) - 1.0/specialpointin );
    }
  }
  
  /* calculate lambda (the negative of the minimun value of f~) */
  lambda = scalar2d_get(f_grid, nz-1, 0, 0);
  for(k=0; k<np; k++) {
    for(j=0; j<nt; j++) {
      /*printf("j=%d, k=%d\n", j, k);*/
      new = scalar2d_get(f_grid, nz-1, j, k);
      lambda = (new < lambda ? new : lambda);
    }
  }
  lambda = -lambda;
  printf("lambda=%f\n", lambda);

  /* calculate alpha */
  alpha = 0.5*(lambda - 1.0/specialpointin);
  gsl_vector_set(alphalist, nz-1, alpha);
  printf("alpha=%f\n", alpha);

  /* rescale f~ to calculate f */ 
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
      scalar2d_set(f_grid, nz-1, j, k, (scalar2d_get(f_grid, nz-1, j, k) + lambda)/alpha);
    }
  }
  
}


/************************************************************************************************/
/* For each point (xi_i, theta_j, phi_k) calculate the radial position in spherical coordinates */
/* using the relations                                                                          */
/*     R(xi, theta, phi) = alpha*[xi + A(xi)*F(theta, phi) + B(xi)*G(theta, phi)]               */ 
/* for the kernel, or                                                                           */
/*     R(xi, theta, phi) = alpha*[xi + A(xi)*F(theta, phi) + B(xi)*G(theta, phi)] + beta        */
/* for the shells, or                                                                           */
/*     U(xi, theta, phi) = 1/r = alpha*[xi + A(xi)*F(theta, phi) - 1]                           */
/* for the external domain.  Return the values of R and U in rgrid.                             */
/************************************************************************************************/
void rofxtp(scalar3d *rgrid, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g)
{
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  double alpha;
  double beta;
  double xi;
  double r; /* the radial position (or 1/r for the external domain) */
  
  nz = rgrid->nz;
  nr = rgrid->nr;
  nt = rgrid->nt;
  np = rgrid->np;

   /* kernel */
  alpha = gsl_vector_get(alphalist, 0);
  for(i=0; i<nr; i++) {
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	r = alpha*(xi 
		   + xi*xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f, 0, j, k)
		   + 0.5*xi*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g, 0, j, k));
	scalar3d_set(rgrid, 0, i, j, k, r);
      }
    }
  }
  
  /* shells */
  for(z=1; z<nz-1; z++) {
    
    alpha = gsl_vector_get(alphalist, z);
    beta = gsl_vector_get(betalist, z);
    for(i=0; i<nr; i++) {
      xi = -cos(PI*i/(nr-1));
      for(j=0; j<nt; j++) {
	for(k=0; k<np; k++) {
	  r = alpha*(xi 
		     + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, z, j, k)
		     + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g, z, j, k))
	    + beta;
	  scalar3d_set(rgrid, z, i, j, k, r);
	}
      }
    }
    
  }
  
  /* external domain */
  alpha = gsl_vector_get(alphalist, nz-1);
  for(i=0; i<nr; i++) {
    xi = -cos(PI*i/(nr-1));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	r = alpha*(xi + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, nz-1, j, k) - 1.0);
	scalar3d_set(rgrid, nz-1, i, j, k, r);
      }
    }
  }
  
}


/* /\*****************************\/ */
/* /\*  d^2 R(xi, theta, phi)    *\/ */
/* /\*  --------------------     *\/ */
/* /\*  d xi^2                   *\/ */
/* /\*****************************\/ */
/* void d2rdxi2(scalar3d *op, scalar2d *f, scalar2d *g) */
/* { */
/*   /\* kernel *\/ */
/*   alpha = ; */
/*   for(i=0; i<nr; i++) { */
/*     xi = ; */
/*     for(j=0; j<nt; j++) { */
/*       for(k=0; k<np; k++) { */
/* 	op = alpha*(xi*xi*(36.0-60.0*xi*xi)*scalar2d_get(f, 0, j, k) */
/* 		    + xi*(15.0-30.0*xi*xi)*scalar2d_get(g, 0, j, k)); */
/* 	scalar3d_set(op, 0, i, j, k); */
/*       } */
/*     } */
/*   } */
  
/*   /\* shells *\/ */
/*   for(z=1; z<nz-1; z++) { */
    
/*     alpha = ; */
/*     beta = ; */
/*     for(i=0; i<nr; i++) { */
/*       xi = ; */
/*       for(j=0; j<nt; j++) { */
/* 	for(k=0; k<np; k++) { */
/* 	  op = alpha*(1.5*xi*scalar2d_get(f, z, j, k) */
/* 		      + 1.5*xi*scalar2d_get(g, z, j, k)); */
/* 	  scalar3d_set(op, z, i, j, k); */
/* 	} */
/*       } */
/*     } */
    
/*   } */
  
/*   /\* external domain *\/ */
/*   alpha = ; */
/*   for(i=0; i<nr; i++) { */
/*     xi = ; */
/*     for(j=0; j<nt; j++) { */
/*       for(k=0; k<np; k++) { */
/* 	op = alpha*1.5*xi*scalar2d_get(f, nz-1, j, k); */
/* 	scalar3d_set(op, nz-1, i, j, k); */
/*       } */
/*     } */
/*   } */
  
/* } */


/* /\**************************************************************\/ */
/* /\* J_1 = dR(xi, theta, phi)/dxi                               *\/ */
/* /\**************************************************************\/ */
/* scalar3d jacobian1(scalar3d *jacobian, scalar2d *f, scalar2d *g) */
/* { */
/*   /\* kernel *\/ */
/*   alpha = ; */
/*   for(i=0; i<nr; i++) { */
/* 	xi = ; */
/* 	for(j=0; j<nt; j++) { */
/* 	  for(k=0; k<np; k++) { */
/* 		jacobian = alpha*(1.0  */
/* 						  + 12.0*xi*xi*xi*(1.0-xi*xi)*scalar2d_get(f, 0, j, k) */
/* 						  + 7.5*xi*xi*(1.0-xi*xi)*scalar2d_get(g, 0, j, k)); */
/* 		scalar3d_set(jacobian, 0, i, j, k); */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* shells *\/ */
/*   for(z=1; z<nz-1; z++) { */
	
/* 	alpha = ; */
/* 	for(i=0; i<nr; i++) { */
/* 	  xi = ; */
/* 	  for(j=0; j<nt; j++) { */
/* 		for(k=0; k<np; k++) { */
/* 		  jacobian = alpha*(1.0  */
/* 							+ 0.75*(xi*xi-1.0)*scalar2d_get(f, z, j, k) */
/* 							+ 0.75*(1.0-xi*xi)*scalar2d_get(g, z, j, k)); */
/* 		  scalar3d_set(jacobian, z, i, j, k); */
/* 		} */
/* 	  } */
/* 	} */
	
/*   } */
  
/*   /\* external domain *\/ */
/*   alpha = ; */
/*   for(i=0; i<nr; i++) { */
/* 	xi = ; */
/* 	for(j=0; j<nt; j++) { */
/* 	  for(k=0; k<np; k++) { */
/* 		jacobian = alpha*(1.0 + 0.75*(xi*xi-1.0)*scalar2d_get(f, nz-1, j, k)); */
/* 		scalar3d_set(jacobian, nz-1, i, j, k); */
/* 	  } */
/* 	} */
/*   } */
  
/* } */


/* /\**************************************************************\/ */
/* /\* J_2 = 1/R * dR(xi, theta, phi)/dtheta                      *\/ */
/* /\**************************************************************\/ */
/* scalar3d jacobian2(scalar3d *jacobian, scalar2d *f, scalar2d *g) */
/* { */
/*   /\* go to Fourier space *\/ */
/*   gridtofourier_bound(f_fourier, f); */
/*   gridtofourier_bound(g_fourier, g); */
  
/*   /\* take d/dtheta derivatives of f and g for all zones *\/ */
/*   /\* and return them in the same structure *\/ */
/*   /\* Derivatives turn sin to cos and vice versa.  *\/ */
/*   /\* d/dt[a_j cos(j theta)] = -j a_j sin(j theta) *\/ */
/*   /\* d/dt[a_j sin(j theta)] = +j a_j cos(j theta) *\/ */
/*   for(z=0; z<nz; z++) { */
/* 	for(k=0; k<npc; k++) { */
/* 	  for(j=0; j<nt; j++) { */
/* 		/\* f: *\/ */
/* 		dcosjtdt = -j * boundcoeff_get(f_fourier, z, j, k, REAL); */
/* 		dsinjtdt = j * boundcoeff_get(f_fourier, z, j, k, IMAG); */
/* 		boundcoeff_set(f_fourier, z, j, k, REAL, dsinjtdt); */
/* 		boundcoeff_set(f_fourier, z, j, k, IMAG, dcosjtdt); */
/* 		/\* g: *\/ */
/* 		dcosjtdt = -j * boundcoeff_get(g_fourier, z, j, k, REAL); */
/* 		dsinjtdt = j * boundcoeff_get(g_fourier, z, j, k, IMAG); */
/* 		boundcoeff_set(g_fourier, z, j, k, REAL, dsinjtdt); */
/* 		boundcoeff_set(g_fourier, z, j, k, IMAG, dcosjtdt); */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* return from Fourier space to grid *\/ */
/*   fouriertogrid_bound(dfdtheta, f_fourier); */
/*   fouriertogrid_bound(dgdtheta, g_fourier); */
  
/*   /\* kernel *\/ */
/*   for(i=0; i<nr; i++) { */
/* 	xi = ; */
/* 	for(j=0; j<nt; j++) { */
/* 	  for(k=0; k<np; k++) { */
/* 		numerator = xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(dfdtheta, 0, j, k) */
/* 		  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(dgdtheta, 0, j, k); */
/* 		denominator = 1.0 + xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f, 0, j, k) */
/* 		  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g, 0, j, k); /\* check the 1/2 factor in front of g *\/ */
/* 		scalar3d_set(jacobian, 0, j, k, numerator/denominator); */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* shells *\/ */
/*   for(z=1; z<nz-1; z++) { */
/* 	alpha = ; */
/* 	beta = ; */
/* 	for(i=0; i<nr; i++) { */
/* 	  xi = ; */
/* 	  for(j=0; j<nt; j++) { */
/* 		for(k=0; k<np; k++) { */
/* 		  numerator = 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(dfdtheta, z, j, k) */
/* 			+ 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(dgdtheta, z, j, k); */
/* 		  denominator =  xi  */
/* 			+ 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, z, j, k) */
/* 			+ 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g, z, j, k) */
/* 			+ beta/alpha; */
/* 		  scalar3d_set(jacobian, z, j, k, numerator/denominator); */
/* 		} */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* external domain *\/ */
/*   for(i=0; i<nr; i++) { */
/* 	xi = ; */
/* 	for(j=0; j<nt; j++) { */
/* 	  for(k=0; k<np; k++) { */
/* 		numerator = 0.25*(xi*xi + xi - 1.0)*scalar2d_get(dfdtheta, nz-1, j, k); */
/* 		denominator = 1.0 + 0.25*(xi*xi + xi - 1.0)*scalar2d_get(f, nz-1, j, k); */
/* 		scalar3d_set(jacobian, nz-1, j, k, numerator/denominator); */
/* 	  } */
/* 	} */
/*   } */
    
/* } */


/* /\**************************************************************\/ */
/* /\* J_3 = 1/Rsin(theta) * dR(xi, theta, phi)/dphi              *\/ */
/* /\**************************************************************\/ */
/* scalar3d jacobian3(scalar3d *jacobian, scalar2d *f, scalar2d *g) */
/* { */

/*   /\* go to Fourier space *\/ */
/*   gridtofourier_bound(f_fourier, f); */
/*   gridtofourier_bound(g_fourier, g); */
 
/*   /\* take d/dtheta derivatives of f and g for all zones *\/ */
/*   /\* and return them in the same structure *\/ */
/*   /\* Derivatives turn sin to cos and vice versa.  *\/ */
/*   /\* d/dp[a_k cos(k phi)] = -k a_k sin(k phi) *\/ */
/*   /\* d/dp[a_k sin(k phi)] = +k a_k cos(k phi) *\/ */
/*   for(z=0; z<nz; z++) { */
/* 	for(k=0; k<npc; k++) { */
/* 	  for(j=0; j<nt; j++) { */
/* 		/\* f: *\/ */
/* 		dcosjtdt = -j * boundcoeff_get(f_fourier, z, j, k, REAL); */
/* 		dsinjtdt = j * boundcoeff_get(f_fourier, z, j, k, IMAG); */
/* 		boundcoeff_set(f_fourier, z, j, k, REAL, dsinjtdt); */
/* 		boundcoeff_set(f_fourier, z, j, k, IMAG, dcosjtdt); */
/* 		/\* g: *\/ */
/* 		dcosjtdt = -j * boundcoeff_get(g_fourier, z, j, k, REAL); */
/* 		dsinjtdt = j * boundcoeff_get(g_fourier, z, j, k, IMAG); */
/* 		boundcoeff_set(g_fourier, z, j, k, REAL, dsinjtdt); */
/* 		boundcoeff_set(g_fourier, z, j, k, IMAG, dcosjtdt); */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* return from Fourier space to grid *\/ */
/*   fouriertogrid_bound(dfdtheta, f_fourier); */
/*   fouriertogrid_bound(dgdtheta, g_fourier); */
  
/*   /\* kernel *\/ */
/*   for(i=0; i<nr; i++) { */
/* 	xi = ; */
/* 	for(j=0; j<nt; j++) { */
/* 	  for(k=0; k<np; k++) { */
/* 		numerator = xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(dfdtheta, 0, j, k) */
/* 		  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(dgdtheta, 0, j, k); */
/* 		denominator = 1.0 + xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f, 0, j, k) */
/* 		  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g, 0, j, k); /\* check the 1/2 factor in front of g *\/ */
/* 		scalar3d_set(jacobian, 0, j, k, numerator/denominator); */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* shells *\/ */
/*   for(z=1; z<nz-1; z++) { */
/* 	alpha = ; */
/* 	beta = ; */
/* 	for(i=0; i<nr; i++) { */
/* 	  xi = ; */
/* 	  for(j=0; j<nt; j++) { */
/* 		for(k=0; k<np; k++) { */
/* 		  numerator = 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(dfdtheta, z, j, k) */
/* 			+ 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(dgdtheta, z, j, k); */
/* 		  denominator =  xi  */
/* 			+ 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, z, j, k) */
/* 			+ 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g, z, j, k) */
/* 			+ beta/alpha; */
/* 		  scalar3d_set(jacobian, z, j, k, numerator/denominator); */
/* 		} */
/* 	  } */
/* 	} */
/*   } */
  
/*   /\* external domain *\/ */
/*   for(i=0; i<nr; i++) { */
/* 	xi = ; */
/* 	for(j=0; j<nt; j++) { */
/* 	  for(k=0; k<np; k++) { */
/* 		numerator = 0.25*(xi*xi + xi - 1.0)*scalar2d_get(dfdtheta, nz-1, j, k); */
/* 		denominator = 1.0 + 0.25*(xi*xi + xi - 1.0)*scalar2d_get(f, nz-1, j, k); */
/* 		scalar3d_set(jacobian, nz-1, j, k, numerator/denominator); */
/* 	  } */
/* 	} */
/*   } */

/* } */

/* /\*******************************************************************\/ */
/* /\* d^2/dtheta^2 + cot(theta) d/dtheta + 1/sin^2(theta) d^2/dphi^2  *\/ */
/* /\*******************************************************************\/ */
/* void laplace_ang(coeff *lapf, coeff *f) */
/* { */


/*   /\* add d^2f/dtheta^2 to lapf *\/ */
/*   for(z=0; z<nz; z++) { */
/* 	for(k=0; k<npc; k++) { */
/* 	  for(j=0; j<nt; j++) { */
/* 		for(i=0; i<nr; i++) { */
/* 		/\* real part *\/ */
/* 		temp = coeff_get(lapf, z, i, j, k, REAL) + -j*j*coeff_get(f, z, i, j, k, REAL); */
/* 		coeff_set(lapf, z, i, j, k, REAL, temp); */
/* 		/\* immaginary part *\/ */
/* 		temp = coeff_get(lapf, z, i, j, k, IMAG) + -j*j*coeff_get(f, z, i, j, k, IMAG); */
/* 		coeff_set(lapf, z, i, j, k, IMAG, temp); */
/* 		} */
/* 	  } */
/* 	} */
/*   } */


/* } */




/**********************************************************************/
/* Take grid points for each boundary b_zjk and decompose             */
/* them into coefficients of {cos(j theta), sin(j theta)}e^(i k phi). */
/*                                                                    */
/* WARNING: THE GRID POINTS CORRESPONDING TO THETA = 0 OR PI MUST BE  */
/*          THE SAME FOR ALL PHI.  OTHERWISE THE GRID IS MULTIVALUED  */
/*          AT THE POLES.  THE INVERSE TRANSFORM WILL NOT RETURN THE  */
/*          ORIGINAL DATA IF THIS CONDITION IS NOT MET.               */
/**********************************************************************/
void gridtofourier_bound(bound_coeff *bcoeff, scalar2d *b_zjk)
{
  int nz; /* number of zones */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int j, k;
  double *in1dphi; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying even theta for fixed phi */
  fftw_plan plan_forward;
  fftw_complex *out1dphic; /* fft in phi direction for fixed theta */
  double *out1dcos; /* fct in theta direction for fixed phi (k is even) */
  double *out1dsin; /* fst in theta direction for fixed phi (k is odd) */
  
  nz = b_zjk->nz;
  nt = b_zjk->nt;
  np = b_zjk->np;
  
  if(np%2 == 1)
	printf("np must be even in gridtofourier_bound\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* Set up arrays to hold the data */
  in1dphi = fftw_malloc ( sizeof ( double ) * np );
  in1dcos = fftw_malloc ( sizeof ( double ) * nt );
  in1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  out1dphic = fftw_malloc ( sizeof ( fftw_complex ) * npc );
  out1dcos = fftw_malloc ( sizeof ( double ) * nt );
  out1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  
  /* big loop for each boundary */
  for (z=0; z<nz; z++) {
	
	/*>>>>>>>>>>>>>>>>>>>>>>> PHI DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<*/

	/* transform for phi for fixed theta value (fixed j) */
	for ( j = 0; j < nt; j++ ) {
	  for ( k = 0; k < np; k++ ) {
		in1dphi[k] = scalar2d_get(b_zjk, z, j, k);
	  }
	  plan_forward = fftw_plan_dft_r2c_1d ( np, in1dphi, out1dphic, FFTW_ESTIMATE );
	  fftw_execute ( plan_forward );
	  for ( k = 0; k < npc; k++ ) {
		bound_coeff_set( bcoeff, z, j, k, REAL,
						 out1dphic[k][0]*(2.0-delta(k, 0)-delta(k, npc-1))/np );
		bound_coeff_set( bcoeff, z, j, k, IMAG,
						 out1dphic[k][1]*(-1)*(2.0-delta(k, 0)-delta(k, npc-1))/np );
	  }
	} 
	/*printf( "Oh hi.\n" );
	  print_bound_coeff( bcoeff );*/

	/*>>>>>>>>>>>>>>>>>>>>>> THETA DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<*/

	/* cosine transform for real part of even k. */
	for ( k = 0; k < npc; k += 2 ) {
	  for ( j = 0; j < nt; j++ ) {
		in1dcos[j] = bound_coeff_get( bcoeff, z, j, k, REAL );
	  }
	  plan_forward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_forward );
	  for ( j = 0; j < nt; j++ ) {
		bound_coeff_set( bcoeff, z, j, k, REAL, 
						 out1dcos[j]*(2.0-delta(j, 0)-delta(j, nt-1))/(2*(nt-1)) );
	  }
	} 
	
	/* cosine transform for immaginary part of even k. */
	/* k = 0 and k = npc-1 don't have immaginary parts */
	for ( k = 2; k < npc-1; k += 2 ) { 
	  for ( j = 0; j < nt; j++ ) {
		in1dcos[j] = bound_coeff_get( bcoeff, z, j, k, IMAG );
	  }
	  plan_forward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_forward );
	  for ( j = 0; j < nt; j++ ) {
		bound_coeff_set( bcoeff, z, j, k, IMAG, 
						 out1dcos[j]*(2.0-delta(j, 0)-delta(j, nt-1))/(2*(nt-1)) );
	  }
	} 
	
	/* sine transform for real part of odd k. */
	for ( k = 1; k < npc; k += 2 ) {
	  /* j = 0 and j = nt-1 don't exist for sine transform */
	  for ( j = 1; j < nt-1; j++ ) {
		in1dsin[j-1] = bound_coeff_get( bcoeff, z, j, k, REAL );
	  }
	  plan_forward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_forward );
	  bound_coeff_set(bcoeff, z, 0, k, REAL, 0.0); /* first element is 0 and not included in RODFT00 */
	  for ( j = 1; j < nt-1; j++ ) {
		bound_coeff_set(bcoeff, z, j, k, REAL, out1dsin[j-1]/(nt-1));
	  }
	  bound_coeff_set(bcoeff, z, nt-1, k, REAL, 0.0); /* last element is 0 and not included in RODFT00 */
	} 
	
	/* sine transform for immaginary part of odd k. */	
	/* k = 0 and k = npc-1 don't have immaginary parts */
	for ( k = 1; k < npc-1; k += 2 ) {
	  /* j = 0 and j = nt-1 don't exist for sine transform */
	  for ( j = 1; j < nt-1; j++ ) {
		in1dsin[j-1] = bound_coeff_get( bcoeff, z, j, k, IMAG );
	  }
	  plan_forward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_forward );
	  bound_coeff_set(bcoeff, z, 0, k, IMAG, 0.0); 
	  for ( j = 1; j < nt-1; j++ ) {
		bound_coeff_set(bcoeff, z, j, k, IMAG, out1dsin[j-1]/(nt-1)); 
	  }
	  bound_coeff_set(bcoeff, z, nt-1, k, IMAG, 0.0); 
	} 	
	
  }
  
  /* Delete plan and arrays */
  fftw_destroy_plan ( plan_forward );
  fftw_free ( in1dphi );
  fftw_free ( out1dphic );
  fftw_free ( in1dcos );
  fftw_free ( in1dsin );
  fftw_free ( out1dcos );
  fftw_free ( out1dsin );
}


/****************************************************************************/
/* Take coefficients of {cos(j theta), sin(j theta)}e^(i k phi)              */
/* and return radial distance to boundary for each (theta, phi).            */
/****************************************************************************/
void fouriertogrid_bound(scalar2d *b_zjk, bound_coeff *bcoeff)
{
  int nz; /* number of zones */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int j, k;
  fftw_complex *in1dphic; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying even theta for fixed phi */
  fftw_plan plan_backward;
  double *out1dphi; /* inverse fft in phi direction for fixed theta */
  double *out1dcos; /* fct in theta direction for fixed phi (k is even) */
  double *out1dsin; /* fst in theta direction for fixed phi (k is odd) */

  nz = b_zjk->nz;
  nt = b_zjk->nt;
  np = b_zjk->np;
  
  if(np%2 == 1)
	printf("np must be even in fouriertogrid_bound\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* Set up arrays to hold the data */
  in1dcos = fftw_malloc ( sizeof ( double ) * nt );
  in1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) ); 
  out1dcos = fftw_malloc ( sizeof ( double ) * nt );
  out1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  in1dphic = fftw_malloc ( sizeof ( fftw_complex ) * npc );
  out1dphi = fftw_malloc ( sizeof ( double ) * np );
 
  
  /* big loop for each boundary */
  for (z=0; z<nz; z++) {
	
	/* cosine transform for real part of even k. */
	for ( k = 0; k < npc; k += 2 ) {
	  for ( j = 0; j < nt; j++ ) {
		in1dcos[j] = bound_coeff_get( bcoeff, z, j, k, REAL ) 
		  / (2.0-delta(j, 0)-delta(j, nt-1)); /* Correct for my convention */
		/*printf("%f  ", in1dcos[j]);*/
	  }
	  /*printf("\n");*/
	  plan_backward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_backward );
	  for ( j = 0; j < nt; j++ ) {
		bound_coeff_set( bcoeff, z, j, k, REAL, out1dcos[j] );
		/*printf("%f  ", out1dcos[j]);*/
	  }
	  /*printf("\n");*/
	} 
	
	/* cosine transform for immaginary part of even k. */
	/* k = 0 and k = npc-1 don't have immaginary parts */
	for ( k = 2; k < npc-1; k += 2 ) { 
	  for ( j = 0; j < nt; j++ ) {
		in1dcos[j] = bound_coeff_get( bcoeff, z, j, k, IMAG )
		  / (2.0-delta(j, 0)-delta(j, nt-1));
	  }
	  plan_backward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_backward );
	  for ( j = 0; j < nt; j++ ) {
		bound_coeff_set( bcoeff, z, j, k, IMAG, out1dcos[j] );
	  }
	} 
	
	/* sine transform for real part of odd k. */
	for ( k = 1; k < npc; k += 2 ) {
	  /* j = 0 and j = nt-1 don't exist for sine transform */
	  for ( j = 1; j < nt-1; j++ ) {
		in1dsin[j-1] = bound_coeff_get( bcoeff, z, j, k, REAL ) / 2.0;
	  }
	  plan_backward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_backward );
	  for ( j = 1; j < nt-1; j++ ) {
		bound_coeff_set(bcoeff, z, j, k, REAL, out1dsin[j-1]);
	  }
	} 
	
	/* sine transform for immaginary part of odd k. */	
	/* k = 0 and k = npc-1 don't have immaginary parts */
	for ( k = 1; k < npc-1; k += 2 ) {
	  /* j = 0 and j = nt-1 don't exist for sine transform */
	  for ( j = 1; j < nt-1; j++ ) {
		in1dsin[j-1] = bound_coeff_get( bcoeff, z, j, k, IMAG ) / 2.0;
	  }
	  plan_backward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
	  fftw_execute ( plan_backward );
	  for ( j = 1; j < nt-1; j++ ) {
		bound_coeff_set(bcoeff, z, j, k, IMAG, out1dsin[j-1]); 
	  }
	} 	
	
	/*printf( "Oh hi.\n" );
	  print_bound_coeff( bcoeff );*/
	
	/* transform for phi for fixed theta value (fixed j) */
	for ( j = 0; j < nt; j++ ) {
	  for ( k = 0; k < npc; k++ ) {
		in1dphic[k][0] = bound_coeff_get(bcoeff, z, j, k, REAL)
		  / (2.0-delta(k, 0)-delta(k, npc-1));
		in1dphic[k][1] = bound_coeff_get(bcoeff, z, j, k, IMAG)
		  * (-1) / (2.0-delta(k, 0)-delta(k, npc-1));
	  }
	  plan_backward = fftw_plan_dft_c2r_1d ( np, in1dphic, out1dphi, FFTW_ESTIMATE );
	  fftw_execute ( plan_backward );
	  for ( k = 0; k < np; k++ ) {
		scalar2d_set( b_zjk, z, j, k, out1dphi[k] );
	  }
	} 
	
	
	
  }
  
  /* Delete plan and arrays */
  fftw_destroy_plan ( plan_backward );
  fftw_free ( in1dphic );
  fftw_free ( out1dphi );
  fftw_free ( in1dcos );
  fftw_free ( in1dsin );
  fftw_free ( out1dcos );
  fftw_free ( out1dsin );
}


/**********************************************************************/
/* Take grid points c_zijk for each zone and decompose them into      */
/* coefficients of T_i(xi){cos(j theta), sin(j theta)}e^(i k phi).    */
/*                                                                    */
/* WARNING: THE GRID POINTS CORRESPONDING TO THETA = 0 OR PI MUST BE  */
/*          THE SAME FOR ALL PHI.  OTHERWISE THE GRID IS MULTIVALUED  */
/*          AT THE POLES.  THE INVERSE TRANSFORM WILL NOT RETURN THE  */
/*          ORIGINAL DATA IF THIS CONDITION IS NOT MET.               */
/**********************************************************************/
void gridtofourier(coeff *coeff, scalar3d *c_zijk)
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int i, j, k;
  double *in1dphi; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying even theta for fixed phi */
  double *in1dxi; /* even part of kernel or all parts of other zones */
  double *in1dxiodd; /* odd part of kernel */
  fftw_plan plan_forward;
  fftw_complex *out1dphic; /* fft in phi direction for fixed theta */
  double *out1dcos; /* fct in theta direction for fixed phi (k is even) */
  double *out1dsin; /* fst in theta direction for fixed phi (k is odd) */
  double *out1dxi;
  double *out1dxiodd;
  int ztest;
  
  nz = c_zijk->nz;
  nr = c_zijk->nr;
  nt = c_zijk->nt;
  np = c_zijk->np;
  
  if(np%2 == 1)
	printf("np must be even in gridtofourier\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* Set up arrays to hold the data */
  in1dphi = fftw_malloc ( sizeof ( double ) * np );
  in1dcos = fftw_malloc ( sizeof ( double ) * nt );
  in1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  in1dxi = fftw_malloc ( sizeof ( double ) * nr );
  in1dxiodd = fftw_malloc ( sizeof ( double ) * (nr-1) );
  out1dphic = fftw_malloc ( sizeof ( fftw_complex ) * npc );
  out1dcos = fftw_malloc ( sizeof ( double ) * nt );
  out1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  out1dxi = fftw_malloc ( sizeof ( double ) * nr );
  out1dxiodd = fftw_malloc ( sizeof ( double ) * (nr-1) );
  
  /* big loop to find coefficients of basis functions for each zone */
  for (z=0; z<nz; z++) {
	
	/*>>>>>>>>>>>>>>>>>>>>>>> PHI DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<*/

	/* transform for phi for fixed xi and theta */
	for ( i = 0; i < nr; i++ ) {
	  for ( j = 0; j < nt; j++ ) {
		for ( k = 0; k < np; k++ ) {
		  in1dphi[k] = scalar3d_get(c_zijk, z, i, j, k);
		}
		plan_forward = fftw_plan_dft_r2c_1d ( np, in1dphi, out1dphic, FFTW_ESTIMATE );
		fftw_execute ( plan_forward );
		for ( k = 0; k < npc; k++ ) {
		  coeff_set( coeff, z, i, j, k, REAL,
					 out1dphic[k][0]*(2.0-delta(k, 0)-delta(k, npc-1))/np );
		  coeff_set( coeff, z, i, j, k, IMAG,
					 out1dphic[k][1]*(-1)*(2.0-delta(k, 0)-delta(k, npc-1))/np );
		}
	  } 
	}
	
	/*>>>>>>>>>>>>>>>>>>>>>> THETA DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<*/

	/* cosine transform for real part of even k. */
	for ( i = 0; i < nr; i++ ) {
	  for ( k = 0; k < npc; k += 2 ) {
		for ( j = 0; j < nt; j++ ) {
		  in1dcos[j] = coeff_get( coeff, z, i, j, k, REAL );
		}
		plan_forward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_forward );
		for ( j = 0; j < nt; j++ ) {
		  coeff_set( coeff, z, i, j, k, REAL, 
					 out1dcos[j]*(2.0-delta(j, 0)-delta(j, nt-1))/(2*(nt-1)) );
		}
	  } 
	}	  
	
	/* cosine transform for immaginary part of even k. */
	for ( i = 0; i < nr; i++ ) {
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 2; k < npc-1; k += 2 ) { 
		for ( j = 0; j < nt; j++ ) {
		  in1dcos[j] = coeff_get( coeff, z, i, j, k, IMAG );
		}
		plan_forward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_forward );
		for ( j = 0; j < nt; j++ ) {
		  coeff_set( coeff, z, i, j, k, IMAG, 
					 out1dcos[j]*(2.0-delta(j, 0)-delta(j, nt-1))/(2*(nt-1)) );
		}
	  } 
	}
	
	/* sine transform for real part of odd k. */
	for ( i = 0; i < nr; i++ ) {
	  for ( k = 1; k < npc; k += 2 ) {
		/* j = 0 and j = nt-1 don't exist for sine transform */
		for ( j = 1; j < nt-1; j++ ) {
		  in1dsin[j-1] = coeff_get( coeff, z, i, j, k, REAL );
		}
		plan_forward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_forward );
		coeff_set(coeff, z, i, 0, k, REAL, 0.0); /* first element is 0 and not included in RODFT00 */
		for ( j = 1; j < nt-1; j++ ) {
		  coeff_set(coeff, z, i, j, k, REAL, out1dsin[j-1]/(nt-1));
		}
		coeff_set(coeff, z, i, nt-1, k, REAL, 0.0); /* last element is 0 and not included in RODFT00 */
	  } 
	}
	
	/* sine transform for immaginary part of odd k. */	
	for ( i = 0; i < nr; i++ ) {
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k += 2 ) {
		/* j = 0 and j = nt-1 don't exist for sine transform */
		for ( j = 1; j < nt-1; j++ ) {
		  in1dsin[j-1] = coeff_get( coeff, z, i, j, k, IMAG );
		}
		plan_forward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_forward );
		coeff_set(coeff, z, i, 0, k, IMAG, 0.0); 
		for ( j = 1; j < nt-1; j++ ) {
		  coeff_set(coeff, z, i, j, k, IMAG, out1dsin[j-1]/(nt-1)); 
		}
		coeff_set(coeff, z, i, nt-1, k, IMAG, 0.0); 
	  } 	
	}
	
/* 	for(ztest=0; ztest<nz; ztest++) {	   */
/* 	  for(i=0; i<nr; i++) { */
/* 		printf("zone %d, xi %d\n\t", z, i); */
/* 		for(k=0; k<npc; k++) */
/* 		  printf("k=%d                   ", k); */
		
/* 		printf("\n"); */
/* 		for ( j = 0; j < nt; j++ ) { */
		  
/* 		  printf("j=%d:\t", j); */
/* 		  for ( k = 0; k < npc; k++ ) { */
/* 			printf ("(%f, %f)  ", coeff_get(coeff, ztest, i, j, k, REAL), coeff_get(coeff, ztest, i, j, k, IMAG)); */
/* 		  } */
/* 		  printf ( "\n" );	    */
		  
/* 		}  */
/* 	  } */
/* 	} */
/* 	printf ( "\n" ); */
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>> XI DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<<<*/


/* 	for ( j = 0; j < nt; j++ ) { */
/* 	  for ( k = 0; k < npc; k++ ) { */
		
/* 		if(!(k%2)) { /\* even k *\/ */
/* 		  /\*printf("in:   ");*\/ */
/* 		  for ( i = 0; i < nr; i++ ) { */
/* 			in1dxi[i] = coeff_get( coeff, z, i, j, k, REAL ); */
/* 			/\*printf("%f ", in1dxi[i]);*\/ */
/* 		  } */
/* 		  /\*printf("\n");*\/ */
/* 		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE ); */
/* 		  fftw_execute ( plan_forward ); */
/* 		  /\*printf("out:  ");*\/ */
/* 		  for ( i = 0; i < nr; i++ ) { */
/* 			coeff_set( coeff, z, i, j, k, REAL,  */
/* 					   out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) ); */
/* 			/\*printf("%f ", coeff_get( coeff, z, i, j, k, REAL ) );	*\/ */
/* 		  } */
/* 		  /\*printf("\n");*\/ */





/* 		} else { /\* odd k *\/ */
/* 		  basis *= sin(j*theta); */
/* 		} */
		
/* 		if( z==0 && !(j%2) ) { /\* even j, zone 0 *\/ */
/* 		  basis *= T_n(2*i, xi); */
/* 		} else if( z==0 && j%2 ) { /\* odd j, zone 0 *\/ */
/* 		  basis *= T_n(2*i+1, xi); */
/* 		} else { /\* other zones *\/ */
/* 		  basis *= T_n(i, xi); */
/* 		} */
		
/* 		sum += basis; */
		
/* 	  } */
/* 	} */
/*   } */



	/* transform for xi for fixed theta and phi */
	if(z == 0) { /* kernel */

	  /* even j ==> even Chebyshev series */
	  /* real part of coefficient */
	  for ( k = 0; k < npc; k++ ) {
		for ( j = 0; j < nt; j += 2 ) {
		  /*printf("in:   ");*/
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, REAL );
			/*printf("%f ", in1dxi[i]);*/
		  }
		  /*printf("\n");*/
		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  /*printf("out:  ");*/
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, REAL, 
					   out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
			/*printf("%f ", coeff_get( coeff, z, i, j, k, REAL ) );	*/
		  }
		  /*printf("\n");*/
		}
	  }
	  /* even j ==> even Chebyshev series */
	  /* immaginary part of coefficient */
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k++ ) {
		for ( j = 0; j < nt; j += 2 ) {
		  /*printf("in:   ");*/
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, IMAG );	
			/*printf("%f ", in1dxi[i]);*/
		  }
		  /*printf("\n");*/
		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  /*printf("out:  ");*/
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG, 
					   out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );	
			/*printf("%f ", coeff_get( coeff, z, i, j, k, IMAG ) );	*/
		  }  
		  /*printf("\n");*/
		}
	  } 
	  /* odd j ==> odd Chebyshev series */
	  /* real part of coefficient */
	  for ( k = 0; k < npc; k++ ) {
		for ( j = 1; j < nt-1; j += 2 ) {
		  /* drop 1st point (it's zero) and reverse other points */		
		  for ( i = 1; i < nr; i++ ) {
			in1dxiodd[i-1] = coeff_get( coeff, z, nr-i, j, k, REAL );
		  }
		  /* do type-3 DCT */
		  plan_forward = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT01, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  for ( i = 0; i < nr-1; i++ ) {
			coeff_set( coeff, z, i, j, k, REAL, 
					   out1dxiodd[i]/(nr-1) );
		  }
		  coeff_set( coeff, z, nr-1, j, k, REAL, 0.0 ); /* set coefficient for T_{nr-1}(xi) = 0 */
		}
	  }
	  /* odd j ==> odd Chebyshev series */
	  /* immaginary part of coefficient */
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k++ ) {
		for ( j = 1; j < nt-1; j += 2 ) {
		  /* drop 1st point (it's zero) and reverse other points */
		  for ( i = 1; i < nr; i++ ) {
			in1dxiodd[i-1] = coeff_get( coeff, z, nr-i, j, k, IMAG );
		  }
		  /* do type-3 DCT */
		  plan_forward = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT01, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  for ( i = 0; i < nr-1; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG, 
					   out1dxiodd[i]/(nr-1) );
		  }  
		  coeff_set( coeff, z, nr-1, j, k, IMAG, 0.0 ); /* set coefficient for T_{nr-1}(xi) = 0 */
		}
	  } 
	  
	} else { /* other zones */
	  
	  /* real part of coefficient */
	  for ( k = 0; k < npc; k++ ) {
		for ( j = 0; j < nt; j++ ) {
		  /*printf("in:   ");*/
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, REAL );
			/*printf("%f ", in1dxi[i]);*/
		  }
		  /*printf("\n");*/
		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  /*printf("out:  ");*/
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, REAL,
					   out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
			/*printf("%f ", coeff_get( coeff, z, i, j, k, REAL ) );*/
		  }
		  /*printf("\n"); */
		}
	  }
	  
	  /* immaginary part of coefficient */
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k++ ) {
		for ( j = 0; j < nt; j++ ) {
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, IMAG );
		  }
		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG,
					   out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
		  }
		  
		}
	  }
	  
	} /* end of radial decomposition */
	
  }
  
  /* Delete plan and arrays */
  fftw_destroy_plan ( plan_forward );
  fftw_free ( in1dphi );
  fftw_free ( out1dphic );
  fftw_free ( in1dcos );
  fftw_free ( in1dsin );
  fftw_free ( out1dcos );
  fftw_free ( out1dsin );
  fftw_free ( in1dxi );
  fftw_free ( out1dxi );
  fftw_free ( in1dxiodd );
  fftw_free ( out1dxiodd );
}


/****************************************************************************/
/* Take coefficients of T_i(xi){cos(j theta), sin(j theta)}e^(i k phi)      */
/* and values on grid (xi_i, theta_j, phi_k).                               */
/****************************************************************************/
void fouriertogrid(scalar3d *c_zijk, coeff *coeff)
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int i, j, k;
  fftw_complex *in1dphic; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying even theta for fixed phi */
  double *in1dxi; /* even part of kernel or all parts of other zones */
  double *in1dxiodd; /* odd part of kernel */
  fftw_plan plan_backward;
  double *out1dphi; /* fft in phi direction for fixed theta */
  double *out1dcos; /* fct in theta direction for fixed phi (k is even) */
  double *out1dsin; /* fst in theta direction for fixed phi (k is odd) */
  double *out1dxi;
  double *out1dxiodd;
  
  nz = c_zijk->nz;
  nr = c_zijk->nr;
  nt = c_zijk->nt;
  np = c_zijk->np;
  
  if(np%2 == 1)
	printf("np must be even in gridtofourier\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* Set up arrays to hold the data */
  in1dphic = fftw_malloc ( sizeof ( fftw_complex ) * npc );
  in1dcos = fftw_malloc ( sizeof ( double ) * nt );
  in1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  in1dxi = fftw_malloc ( sizeof ( double ) * nr );
  in1dxiodd = fftw_malloc ( sizeof ( double ) * (nr-1) );
  out1dphi = fftw_malloc ( sizeof ( double ) * np );
  out1dcos = fftw_malloc ( sizeof ( double ) * nt );
  out1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  out1dxi = fftw_malloc ( sizeof ( double ) * nr );
  out1dxiodd = fftw_malloc ( sizeof ( double ) * (nr-1) );
  
  /* big loop to find coefficients of basis functions for each zone */
  for (z=0; z<nz; z++) {

	/*>>>>>>>>>>>>>>>>>>>>>> INVERSE XI TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/

	/* transform for xi for fixed theta and phi */
	if(z == 0) { /* kernel */

	  /* even j ==> even Chebyshev series */
	  /* real part of coefficient */
	  for ( k = 0; k < npc; k++ ) {
		for ( j = 0; j < nt; j += 2 ) {
		  /*printf("in:   ");*/
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, REAL )
			  / (neg1toi(i)*(2-delta(i, 0)-delta(i, nr-1)));
			/*printf("%f ", in1dxi[i]);*/
		  }
		  /*printf("\n");*/
		  plan_backward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_backward );
		  /*printf("out:  ");*/
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, REAL, out1dxi[i] );
			/*printf("%f ", coeff_get( coeff, z, i, j, k, REAL ) );	*/
		  }
		  /*printf("\n");*/
		}
	  }
	  /* even j ==> even Chebyshev series */
	  /* immaginary part of coefficient */
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k++ ) {
		for ( j = 0; j < nt; j += 2 ) {
		  /*printf("in:   ");*/
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, IMAG ) 
			  / (neg1toi(i)*(2-delta(i, 0)-delta(i, nr-1)));	
			/*printf("%f ", in1dxi[i]);*/
		  }
		  /*printf("\n");*/
		  plan_backward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_backward );
		  /*printf("out:  ");*/
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG, out1dxi[i] );	
			/*printf("%f ", coeff_get( coeff, z, i, j, k, IMAG ) );	*/
		  }  
		  /*printf("\n");*/
		}
	  } 
	  /* odd j ==> odd Chebyshev series */
	  /* real part of coefficient */
	  for ( k = 0; k < npc; k++ ) {
		for ( j = 1; j < nt-1; j += 2 ) { 
		  /*printf("in:   ");*/
		  for ( i = 0; i < nr-1; i++ ) {
			in1dxiodd[i] = coeff_get( coeff, z, i, j, k, REAL );
			/*printf("%f ", in1dxiodd[i]);*/
		  }
		  /*printf("\n");*/
		  /* do type-2 DCT (inverse of type-3 DCT) */
		  plan_backward = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT10, FFTW_ESTIMATE );
		  fftw_execute ( plan_backward );
		  /*printf("out:  ");*/
		  /* 1st point is zero */ 
		  coeff_set( coeff, z, 0, j, k, REAL, 0.0 );
		  /* reverse other points */		
		  for ( i = 0; i < nr-1; i++ ) {
			/*printf("%f ", out1dxiodd[i] );*/
			coeff_set( coeff, z, i+1, j, k, REAL, out1dxiodd[nr-2-i] / 2.0 );
			/*printf("%f ", coeff_get( coeff, z, i+1, j, k, REAL ) );*/
		  }
		  /*intf("\n");*/
		}
	  }
	  /* odd j ==> odd Chebyshev series */
	  /* immaginary part of coefficient */
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k++ ) {
		for ( j = 1; j < nt-1; j += 2 ) {
		  for ( i = 0; i < nr-1; i++ ) {
			in1dxiodd[i] = coeff_get( coeff, z, i, j, k, IMAG );
		  }
		  /* do type-2 DCT (inverse of type-3 DCT) */
		  plan_backward = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT10, FFTW_ESTIMATE );
		  fftw_execute ( plan_backward );
		  /* 1st point is zero */ 
		  coeff_set( coeff, z, 0, j, k, IMAG, 0.0 );
		  /* reverse other points */
		  for ( i = 0; i < nr-1; i++ ) {
			coeff_set( coeff, z, i+1, j, k, IMAG, out1dxiodd[nr-2-i] / 2.0 );
		  }  
		}
	  } 
	  
	} else { /* other zones */
	  
	  /* real part of coefficient */
	  for ( k = 0; k < npc; k++ ) {
		for ( j = 0; j < nt; j++ ) {
		  /*printf("in:   ");*/
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, REAL ) 
			  / (neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1)));
			/*printf("%f ", in1dxi[i]);*/
		  }
		  /*printf("\n");*/
		  plan_backward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_backward );
		  /*printf("out:  ");*/
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, REAL, out1dxi[i] );
			/*printf("%f ", coeff_get( coeff, z, i, j, k, REAL ) );*/
		  }
		  /*printf("\n"); */
		}
	  }
	  
	  /* immaginary part of coefficient */
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k++ ) {
		for ( j = 0; j < nt; j++ ) {
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, IMAG )
			  / (neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1)));
		  }
		  plan_backward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_backward );
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG, out1dxi[i] );
		  }
		}
	  }
	  
	} 

	/*>>>>>>>>>>>>>>>>>>>>>> INVERSE THETA TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/

	/* cosine transform for real part of even k. */
	for ( i = 0; i < nr; i++ ) {
	  for ( k = 0; k < npc; k += 2 ) {
		for ( j = 0; j < nt; j++ ) {
		  in1dcos[j] = coeff_get( coeff, z, i, j, k, REAL ) 
			/ (2.0-delta(j, 0)-delta(j, nt-1)); /* Correct for my convention */
		  /*printf("%f  ", in1dcos[j]);*/
		}
		/*printf("\n");*/
		plan_backward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_backward );
		for ( j = 0; j < nt; j++ ) {
		  coeff_set( coeff, z, i, j, k, REAL, out1dcos[j] );
		  /*printf("%f  ", out1dcos[j]);*/
		}
		/*printf("\n");*/
	  } 
	}
	/* cosine transform for immaginary part of even k. */
	for ( i = 0; i < nr; i++ ) {
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 2; k < npc-1; k += 2 ) { 
		for ( j = 0; j < nt; j++ ) {
		  in1dcos[j] = coeff_get( coeff, z, i, j, k, IMAG )
			/ (2.0-delta(j, 0)-delta(j, nt-1));
		}
		plan_backward = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_backward );
		for ( j = 0; j < nt; j++ ) {
		  coeff_set( coeff, z, i, j, k, IMAG, out1dcos[j] );
		}
	  } 
	} 
	/* sine transform for real part of odd k. */
	for ( i = 0; i < nr; i++ ) {
	  for ( k = 1; k < npc; k += 2 ) {
		/* j = 0 and j = nt-1 don't exist for sine transform */
		for ( j = 1; j < nt-1; j++ ) {
		  in1dsin[j-1] = coeff_get( coeff, z, i, j, k, REAL ) / 2.0;
		}
		plan_backward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_backward );
		for ( j = 1; j < nt-1; j++ ) {
		  coeff_set(coeff, z, i, j, k, REAL, out1dsin[j-1]);
		}
	  } 
	}  
	/* sine transform for immaginary part of odd k. */	
	for ( i = 0; i < nr; i++ ) {
	  /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( k = 1; k < npc-1; k += 2 ) {
		/* j = 0 and j = nt-1 don't exist for sine transform */
		for ( j = 1; j < nt-1; j++ ) {
		  in1dsin[j-1] = coeff_get( coeff, z, i, j, k, IMAG ) / 2.0;
		}
		plan_backward = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
		fftw_execute ( plan_backward );
		for ( j = 1; j < nt-1; j++ ) {
		  coeff_set(coeff, z, i, j, k, IMAG, out1dsin[j-1]); 
		}
	  } 	
	}
 
	
	/*>>>>>>>>>>>>>>>>>>>>>> INVERSE PHI TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/
 
	/* transform for phi for fixed theta value (fixed j) */
	for ( i = 0; i < nr; i++ ) {
	  for ( j = 0; j < nt; j++ ) {
		for ( k = 0; k < npc; k++ ) {
		  in1dphic[k][0] = coeff_get(coeff, z, i, j, k, REAL)
			/ (2.0-delta(k, 0)-delta(k, npc-1));
		  in1dphic[k][1] = coeff_get(coeff, z, i, j, k, IMAG)
			* (-1) / (2.0-delta(k, 0)-delta(k, npc-1));
		}
		plan_backward = fftw_plan_dft_c2r_1d ( np, in1dphic, out1dphi, FFTW_ESTIMATE );
		fftw_execute ( plan_backward );
		for ( k = 0; k < np; k++ ) {
		  scalar3d_set( c_zijk, z, i, j, k, out1dphi[k] );
		}
	  } 
	}
	
  } 

  /* Delete plan and arrays */
  fftw_destroy_plan ( plan_backward );
  fftw_free ( in1dphic );
  fftw_free ( out1dphi );
  fftw_free ( in1dcos );
  fftw_free ( in1dsin );
  fftw_free ( out1dcos );
  fftw_free ( out1dsin );
  fftw_free ( in1dxi );
  fftw_free ( out1dxi );
  fftw_free ( in1dxiodd );
  fftw_free ( out1dxiodd );
}






/********************************************************************/
/* Set grid points.                                                 */
/* Currently a gross waste of memory. Change the structure form.    */
/* Although for the physical grid, r depends on theta and phi.      */
/********************************************************************/
void makegrid3d(grid3d *x_zijk)
{
  int nz; /* number of zones */
  int nr; /* number of radial points per zone */ 
  int nt; /* number of theta points */
  int np; /* number of phi points */ 
  int z; /* current zone */
  int i; /* current radial index */
  int j; /* current theta index */
  int k; /* current phi index */
  
  nz = x_zijk->nz;
  nr = x_zijk->nr;
  nt = x_zijk->nt;
  np = x_zijk->np;
  
  /* set each point */
  for(z=0; z<nz; z++){
	for(i=0; i<nr; i++){
	  for(j=0; j<nt; j++){
		for(k=0; k<np; k++){
		  if(z==0)
			grid3d_set(x_zijk, z, i, j, k, 
					   sin(PI*i/(2*(nr-1))), 
					   PI*j/(nt-1),
					   2*PI*k/np);
		  else
			grid3d_set(x_zijk, z, i, j, k, 
					   -cos(PI*i/(nr-1)),
					   PI*j/(nt-1),
					   2*PI*k/np);
		}
	  }
	}
  }
}


/*******************************************************************/
/* Take an analytical function for the boundary r=S(z, theta, phi) */
/* and evaluate it at the collocation points.                      */
/* Return the value r at those points in b_zjk.                    */
/*******************************************************************/
void boundtogrid(scalar2d *b_zjk, grid3d *points, double (*bound)(int, double, double))
{
  int nz; /* number of zones */
  int nt; /* number of theta points */
  int np; /* number of phi points */ 
  int z; /* current zone */
  int j; /* current theta index */
  int k; /* current phi index */
  double theta;
  double phi;
  
  nz = b_zjk->nz;
  nt = b_zjk->nt;
  np = b_zjk->np;
  
  for(z=0; z<nz; z++){
	for(j=0; j<nt; j++){
	  for(k=0; k<np; k++){
		/* i=1 is arbitrary. don't care about radial part. */
		theta = grid3d_get(points, z, 1, j, k, 1); 
		phi = grid3d_get(points, z, 1, j, k, 2); 
		/*printf("(%d, %f, %f) ", z, theta, phi);*/
		scalar2d_set(b_zjk, z, j, k, bound(z, theta, phi));
		/*printf("%f\t", scalar2d_get(b_zjk, z, j, k));*/
		/*printf("%f\t", bound(z, theta, phi));*/
	  }
	  /*printf("\n");*/
	}
  }
}


/*********************************************************************/
/* Take an analytical function f(z, r, theta, phi)                   */
/* and evaluate it at the collocation points (xi_i, theta_j, phi_k). */
/* Return the value of f at those points in f_zijk.                  */
/*********************************************************************/
void ftogrid(scalar3d *f_zijk, scalar2d *b_zjk, grid3d *points, double (*func)(int, double, double, double))
{
  int nz; /* number of zones */
  int nr; /* number of xi points in each zone */
  int nt; /* number of theta points */
  int np; /* number of phi points */ 
  int z; /* current zone */
  int i; /* current xi index */
  int j; /* current theta index */
  int k; /* current phi index */
  double xi;
  double theta;
  double phi;
  double rin;
  double rout;
  double r;
  
  nz = points->nz;
  nr = points->nr; 
  nt = points->nt;
  np = points->np;
  
  for(j=0; j<nt; j++){
	for(k=0; k<np; k++){
	  theta = grid3d_get(points, 0, 0, j, k, 1); /* z, i not set yet */
	  phi = grid3d_get(points, 0, 0, j, k, 2); 
	  /* theta, phi now fixed. Evaluate along ray from r=0 to infinity: */
	  
	  /* evaluate in kernel */
	  z=0;
	  for(i=0; i<nr; i++){
		xi = grid3d_get(points, z, i, j, k, 0);
		rout = scalar2d_get(b_zjk, z, j, k);
		r = rout*xi;
		scalar3d_set(f_zijk, z, i, j, k, func(z, r, theta, phi));
	  }
	  
	  /* evaluate in shells */
	  for(z=1; z<nz-1; z++){
		for(i=0; i<nr; i++){
		  xi = grid3d_get(points, z, i, j, k, 0);
		  rin = scalar2d_get(b_zjk, z-1, j, k);
		  rout = scalar2d_get(b_zjk, z, j, k);
		  r = 0.5*(rout - rin)*xi + 0.5*(rout + rin);
		  scalar3d_set(f_zijk, z, i, j, k, func(z, r, theta, phi));
		}
	  }
  
	  /* evaluate in external domain */
	  z=nz-1;
	  /* i=nr-1 corresponds to r=infinity. Don't evaluate it. */
	  for(i=0; i<nr-1; i++){
		xi = grid3d_get(points, z, i, j, k, 0);
		rin = scalar2d_get(b_zjk, z-1, j, k);
		r = 2.0*rin/(1-xi);
		scalar3d_set(f_zijk, z, i, j, k, func(z, r, theta, phi));
	  }
	  /* explicitly set f=0 at r=infinity (i=nr-1) */
	  scalar3d_set(f_zijk, z, nr-1, j, k, 0);
	}
  }
}
