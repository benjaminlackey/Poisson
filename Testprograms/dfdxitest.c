/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/poisson.h dfdxitest.c */


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
/*double fieldphysical(int z, double r, double theta, double phi);*/
double T_n(int n, double x);
void fouriertogrid_xideriv(scalar3d *c_zijk, coeff *coeff, int xishift, int thetashift, int phishift);

int main (void)
{

 
  /*clock_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12; 
    float ratio;*/ 

  /*FILE *fpgrid;
    FILE *fpfunction;*/

  int z, i, j, k;
  double xi_i, theta_j, phi_k;
  int nz = 3;
  int nr = 5; /* must be odd? */
  int nt = 7;
  int np = 6; /* must be even */
  int npc;
  scalar2d *b_zjk;
  bound_coeff *bcoeff_zjk;
  scalar3d *c_zijk;
  coeff *coeff_zijk;
  coeff *dfdxi_coeff;
  scalar3d *dfdxi_grid;
  
  
  b_zjk = scalar2d_alloc(nz-1, nt, np);
  bcoeff_zjk = bound_coeff_alloc(nz-1, nt, np);
  
  c_zijk = scalar3d_alloc(nz, nr, nt, np);
  coeff_zijk = coeff_alloc(nz, nr, nt, np);

  dfdxi_grid = scalar3d_alloc(nz, nr, nt, np);
  dfdxi_coeff = coeff_alloc(nz, nr, nt, np);
  
  npc = ( np / 2 ) + 1;
  
  /* Assign the field data. */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if(z==0) {
	    xi_i = sin(PI*i/(2*(nr-1)));
	  } else {
	    xi_i = -cos(PI*i/(nr-1));
	  }
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  scalar3d_set(c_zijk, z, i, j, k, field(z, nr, nt, np, xi_i, theta_j, phi_k));
	}
      }
    }
  }
  
  /*print_scalar3d(c_zijk);*/
 
  gridtofourier(coeff_zijk, c_zijk);
  printf("coefficients before derivative\n");
  print_coeff(coeff_zijk);

  dfdxi(dfdxi_coeff, coeff_zijk);
  printf("coefficients after derivative\n");
  print_coeff(dfdxi_coeff);
  
  fouriertogrid(dfdxi_grid, dfdxi_coeff);
  /*fouriertogrid_xideriv(dfdxi_grid, dfdxi_coeff, 1, 0, 0);*/
 
  gridtofourier(coeff_zijk, c_zijk);

  return 0;
}


/*******************************************************************/
/* return function value at the grid point located at (theta, phi) */
/*******************************************************************/
/* double boundary(int nt, int np, int z, double theta, double phi) */
/* { */
/*   if(z==0) */
/*     return 1.0*(1.0 + 0.1*sin(theta)*(cos(phi)+sin(phi)) + 0.1*(1-cos(2*theta))*(cos(2*phi)+sin(2*phi))); */
/*   else */
/*     return 5.0; */
/* } */
double boundary(int z, int nt, int np, double theta, double phi)
{
  int j, k;
  double sum;
  int npc = np/2 + 1;

  sum = 1;
  for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < npc; k+=2 ) {
	  sum += /*((j+1)*10+(k+1))*/1.1111111*cos(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
	}
	for ( k = 1; k < npc; k+=2 ) {
	  sum += /*((j+1)*10+(k+1))*/1.1111111*sin(j*theta)*(cos(k*phi) + (1-delta(k, 0)-delta(k, npc-1))*sin(k*phi));
	}
  }
  
  /*sum = 0;
  for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np/2+1; k+=2 ) {
	  sum += (j+1)*cos(j*theta)*(cos(k*phi) + sin(k*phi));
	}
	for ( k = 1; k < np/2+1; k+=2 ) {
	  sum += (j+1)*sin(j*theta)*(cos(k*phi) + sin(k*phi));
	}
	}*/
  
  /*sum = 0;
	for ( j = 0; j < nt; j++ ) {
	sum += (j+1)*cos(j*theta);
	}*/
  
  return sum;
}


/*********************************************************/
/* Some function of position in spherical coordinates.   */
/*********************************************************/
/* double fieldphysical(int z, double r, double theta, double phi) */
/* { */
/*   return cos(r); */
/* } */


/******************************************************/
/* Evaluate a function at the point (xi, theta, phi). */
/******************************************************/
double field(int z, int nr, int nt, int np, double xi, double theta, double phi)
{
  int i;
  double sum;
  
  if(z==0) {
    sum = 0.0;
    for(i=0; i<nr; i++){
      sum += T_n(2*i, xi);
    }
    for(i=0; i<nr-1; i++){
      sum += T_n(2*i+1, xi)*(cos(theta) + sin(theta)*(cos(phi)+sin(phi)));
    }
    return sum;
  } else {
    sum = 0.0;
    for(i=0; i<nr; i++){
      sum += T_n(i, xi);
    }
    return sum;
  }
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


/* void coeff_to_even_or_all_cheb(z, j, k, nr, in1dxi, out1dxi, plan_backward_xi) */
/* { */
/*   int i; */
/*   for ( i = 0; i < nr; i++ ) { */
/*     in1dxi[i] = coeff_get( coeff, z, i, j, k, REAL ) */
/*       / (neg1toi(i)*(2-delta(i, 0)-delta(i, nr-1))); */
/*   } */
/*   /\*printf("\n");*\/ */
/*   plan_backward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE ); */
/*   fftw_execute ( plan_backward ); */
/*   /\*printf("out:  ");*\/ */
/*   for ( i = 0; i < nr; i++ ) { */
/*     coeff_set( coeff, z, i, j, k, REAL, out1dxi[i] ); */
/*     /\*printf("%f ", coeff_get( coeff, z, i, j, k, REAL ) );	*\/ */
/*   } */
/*   /\*printf("\n");*\/ */
/* } */
/* } */


/****************************************************************************/
/* Take coefficients of T_i(xi){cos(j theta), sin(j theta)}e^(i k phi)      */
/* and values on grid (xi_i, theta_j, phi_k).                               */
/****************************************************************************/
void fouriertogrid_xideriv(scalar3d *c_zijk, coeff *coeff, int xishift, int thetashift, int phishift)
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int imag;
  int i, j, k;
  fftw_complex *in1dphic; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying odd theta for fixed phi */
  double *in1dxi; /* even part of kernel or all parts of other zones */
  double *in1dxiodd; /* odd part of kernel */
  fftw_plan plan_backward_xi;
  fftw_plan plan_backward_xiodd;
  fftw_plan plan_backward_cos;
  fftw_plan plan_backward_sin;
  fftw_plan plan_backward_phi;
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
  
  /* Allocate memory for 1 dimensional transform arrays. */
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
  
  /* Create FFTW plans.  The algorithm will be stored in the fftw_plan structure.         */
  /* This structure contains pointers to the in and out arrays right?                     */
  /* The arrays can be changed later (must stay same size) but the plan will be the same. */
  /* The plan usually overwrites the data, so set the data after making the plan.         */
  /* It might be faster to figure out how to use FFTW wisdom to speed up planning.        */
  plan_backward_xi = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
  plan_backward_xiodd = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT10,  FFTW_ESTIMATE); /* do type-2 DCT (inverse of type-3 DCT) */
  plan_backward_cos = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
  plan_backward_sin = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
  plan_backward_phi = fftw_plan_dft_c2r_1d ( np, in1dphic, out1dphi, FFTW_ESTIMATE );
  
  
  /* big loop to find coefficients of basis functions for each zone */
  for (z=0; z<nz; z++) {
    
    /*>>>>>>>>>>>>>>>>>>>>>> INVERSE XI TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/
    /*           transform for xi for fixed theta and phi                 */
    
    if(z == 0) { /* xi transform in kernel is either even or odd */
      
      /* even Chebyshev series if j is even or */
      /* j is odd and coeff was acted on by an odd number of d/dxi derivitives */
      for(imag=0; imag<=1; imag++) {
	for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */ 
	  for(j = xishift; j < nt-xishift; j += 2) {
	    for ( i = 0; i < nr; i++ ) {
	      in1dxi[i] = coeff_get( coeff, z, i, j, k, imag ) / (neg1toi(i)*(2-delta(i, 0)-delta(i, nr-1)));
	    }
	    fftw_execute ( plan_backward_xi );
	    for ( i = 0; i < nr; i++ ) {
	      coeff_set( coeff, z, i, j, k, imag, out1dxi[i] );
	    }
	  }
	}
      }
      
      /* odd Chebyshev series if j is odd or */
      /* j is even and coeff was acted on by an odd number of d/dxi derivitives */
      for(imag=0; imag<=1; imag++) {
	for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */ 
	  for(j = 1-xishift; j < nt-1+xishift; j += 2) {
	    for ( i = 0; i < nr-1; i++ ) {
	      in1dxiodd[i] = coeff_get( coeff, z, i, j, k, imag );
	    }
	    /* do type-2 DCT (inverse of type-3 DCT) */
	    fftw_execute ( plan_backward_xiodd );
	    /* 1st point is zero */ 
	    coeff_set( coeff, z, 0, j, k, imag, 0.0 );
	    /* reverse other points */		
	    for ( i = 0; i < nr-1; i++ ) {	
	      coeff_set( coeff, z, i+1, j, k, imag, out1dxiodd[nr-2-i] / 2.0 );	
	    }
	  }
	}
      }
	  
    } else { /* other zones */
      
      for(imag=0; imag<=1; imag++) {
	for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */ 
	  for ( j = 0; j < nt; j++ ) {
	    for ( i = 0; i < nr; i++ ) {
	      in1dxi[i] = coeff_get( coeff, z, i, j, k, imag ) / (neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1)));
	    }
	    fftw_execute ( plan_backward_xi );
	    for ( i = 0; i < nr; i++ ) {
	      coeff_set( coeff, z, i, j, k, imag, out1dxi[i] );
	    }
	  }
	}
      }
      
    } 
    
    /*>>>>>>>>>>>>>>>>>>>>>> INVERSE THETA TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /* Cosine series if k is even. */
    for(imag=0; imag<=1; imag++) {
      for ( i = 0; i < nr; i++ ) {
	for ( k = 2*imag; k < npc-imag; k += 2 ) { /* start k at 0 for real part and 2 for imag part */
	  for ( j = 0; j < nt; j++ ) {
	    in1dcos[j] = coeff_get( coeff, z, i, j, k, imag ) / (2.0-delta(j, 0)-delta(j, nt-1)); 
	  }
	  fftw_execute ( plan_backward_cos );
	  for ( j = 0; j < nt; j++ ) {
	    coeff_set( coeff, z, i, j, k, imag, out1dcos[j] );
	  }
	} 
      }
    }
    
    /* Sine series if k is odd. */
    for(imag=0; imag<=1; imag++) {
      for ( i = 0; i < nr; i++ ) {
	for ( k = 1; k < npc-imag; k += 2 ) {
	  /* j = 0 and j = nt-1 don't exist for sine transform */
	  for ( j = 1; j < nt-1; j++ ) {
	    in1dsin[j-1] = coeff_get( coeff, z, i, j, k, imag ) / 2.0;
	  }
	  fftw_execute ( plan_backward_sin );
	  for ( j = 1; j < nt-1; j++ ) {
	    coeff_set(coeff, z, i, j, k, imag, out1dsin[j-1]);
	  }
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
	fftw_execute ( plan_backward_phi );
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set( c_zijk, z, i, j, k, out1dphi[k] );
	}
      } 
    }
    
  } /* end of z loop */

  /* Delete plan and arrays */
  fftw_destroy_plan ( plan_backward_xi );
  fftw_destroy_plan ( plan_backward_xiodd );
  fftw_destroy_plan ( plan_backward_cos );
  fftw_destroy_plan ( plan_backward_sin );
  fftw_destroy_plan ( plan_backward_phi );
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
