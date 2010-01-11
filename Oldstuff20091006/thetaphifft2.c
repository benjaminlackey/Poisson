/* To compile type: gcc -lm -lfftw3 -lgsl -lgslcblas thetaphifft2.c print.c coefficients.c poisson.h/*


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

double boundary(int nt, int np, double theta, double phi);
void gridtofourier_bound(bound_coeff *bcoeff_zjk, scalar2d *b_zjk);
void fouriertogrid_bound(scalar2d *b_zjk, bound_coeff *bcoeff);
void gridtofourier(coeff *coeff, scalar3d *c_zijk);
double field(int z, int nr, int nt, int np, double xi, double theta, double phi);
double T_n(int n, double x);

int main (void)
{
  int z, i, j, k;
  double xi_i, theta_j, phi_k;
  int nz = 2;
  int nr = 5; /* must be odd? */
  int nt = 7;
  int np = 6; /* must be even */
  int npc;
  scalar2d *b_zjk;
  bound_coeff *bcoeff_zjk;
  scalar3d *c_zijk;
  coeff *coeff_zijk;

  b_zjk = scalar2d_alloc(nz-1, nt, np);
  bcoeff_zjk = bound_coeff_alloc(nz-1, nt, np);

  c_zijk = scalar3d_alloc(nz, nr, nt, np);
  coeff_zijk = coeff_alloc(nz, nr, nt, np);

  npc = ( np / 2 ) + 1;
  
  /* Assign the data. */
  for ( z = 0; z < nz-1; z++ ) {
	for ( k = 0; k < np; k++ ) {
	  phi_k = 2*PI*k/np;
	  for ( j = 0; j < nt; j++ ) {
		theta_j = PI*j/(nt-1);
		scalar2d_set(b_zjk, z, j, k, boundary(nt, np, theta_j, phi_k));
	  }
	}
  }
  
/*   /\* Assign the boundary data version 2. *\/ */
/*   for ( z = 0; z < nz-1; z++ ) { */
/* 	for ( k = 0; k < np; k++ ) { */
/* 	  scalar2d_set(b_zjk, z, 0, k, (z+1)*100); */
/* 	  for ( j = 1; j < nt-1; j++ ) { */
/* 		scalar2d_set(b_zjk, z, j, k, (z+1)*100 + (j+1)*10 + (k+1)); */
/* 	  } */
/* 	  scalar2d_set(b_zjk, z, nt-1, k, (z+1)*100); */
/* 	} */
/*   } */
  
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
  
  /* print start data*/
  printf("  Boundary data.\n");
  
  print_scalar2d(b_zjk);
  
  gridtofourier_bound(bcoeff_zjk, b_zjk);
  
  print_bound_coeff(bcoeff_zjk);
  
  fouriertogrid_bound(b_zjk, bcoeff_zjk);

  print_scalar2d(b_zjk);

  printf("  Field data.\n");

  print_scalar3d(c_zijk);

  gridtofourier(coeff_zijk, c_zijk);

  print_coeff_2(coeff_zijk);

  return 0;
}


double boundary(int nt, int np, double theta, double phi)
{
  int j, k;
  double sum;
  int npc = np/2 + 1;

  sum = 0;
  for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < npc; k+=2 ) {
	  sum += ((j+1)*10+(k+1))*cos(j*theta)*(cos(k*phi) + sin(k*phi));
	}
	for ( k = 1; k < npc-1; k+=2 ) {
	  sum += ((j+1)*10+(k+1))*sin(j*theta)*(cos(k*phi) + sin(k*phi));
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


double field(int z, int nr, int nt, int np, double xi, double theta, double phi)
{
  int i, j, k;
  double fac; /* basis function */
  double sum;
  int npc = np/2 + 1;
  
  sum = 0;
  for ( i = 0; i < nr; i++ ) { 
	for ( j = 0; j < nt; j++ ) {
	  for ( k = 0; k < npc; k++ ) {
		
		fac = cos(k*phi) + sin(k*phi);
		
		if(!(k%2)) { /* even k */
		  fac *= cos(j*theta);
		} else { /* odd k */
		  fac *= sin(j*theta);
		}
		
		if( z==0 && !(j%2) ) { /* even j, zone 0 */
		  fac *= T_n(2*i, xi);
		} else if( z==0 && j%2 ) { /* odd j, zone 0 */
		  fac *= T_n(2*i+1, xi);
		} else { /* other zones */
		  fac *= T_n(i, xi);
		}
		
		sum += 5*fac;
		
	  }
	}
  }
  
/*   sum = 0;   */
/*   if(z==0) { */
/* 	/\* kernel *\/ */
/* 	/\* even k ==> ( cos(j theta) ) *\/ */
/* 	for ( k = 0; k < npc; k+=2 ) { */
/* 	  for ( j = 0; j < nt; j += 2 ) { */
/* 		for ( i = 0; i < nr; i++ ) {  */
/* 		  if  */
/* 		  sum += 7*T_n(2*i, xi)*cos(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 		}  */
/* 	  } */
/* 	} */
/* 	/\* odd k ==> ( sin(j theta) ) *\/ */
/* 	for ( k = 1; k < npc-1; k+=2 ) { */
/* 	  for ( j = 1; j < nt-1; j+=2 ) { */
/* 		for ( i = 1; i < nr-1; i++ ) {  */
/* 		  sum += 7*T_n(2*i+1, xi)*sin(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 		} */
/* 	  } */
/* 	} */
	
/*   } else { */
/* 	/\* other zones *\/ */
/* 	/\* even k ==> even j ( cos(j theta) ) *\/ */
/* 	for ( k = 0; k < npc; k+=2 ) { */
/* 	  for ( j = 0; j < nt; j+=2 ) { */
/* 		for ( i = 0; i < nr; i++ ) {  */
/* 		  sum += 7*T_n(i, xi)*cos(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 		}  */
/* 	  } */
/* 	}	 */
/* 	/\* odd k ==> odd j ( sin(j theta) ) *\/ */
/* 	for ( k = 1; k < npc-1; k+=2 ) { */
/* 	  for ( j = 1; j < nt-1; j+=2 ) { */
/* 		for ( i = 1; i < nr-1; i++ ) {  */
/* 		  sum += 7*T_n(i, xi)*sin(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 		} */
/* 	  } */
/* 	} */
/*   } */
  

/*   sum = 0; */
/*   for ( j = 0; j < nt; j++ ) { */
/* 	for ( k = 0; k < npc; k+=2 ) { */
/* 	  sum += ((j+1)*10+(k+1))*cos(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	} */
/* 	for ( k = 1; k < npc-1; k+=2 ) { */
/* 	  sum += ((j+1)*10+(k+1))*sin(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	}  */
/*   } */
  
  /* sum = 0;  */
/*   for ( i = 0; i < nr; i++ ) { */
/* 	sum += T_n(i, xi); */
/*   } */


  /*   /\* even Chebyshev polynomials *\/ */
/*   sum = 0; */
/*   for ( i = 0; i < nr; i += 2 ) { */
/* 	for ( j = 0; j < nt; j++ ) { */
/* 	  for ( k = 0; k < npc; k+=2 ) { */
/* 		sum += T_n(i, xi)*cos(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	  } */
/* 	  for ( k = 1; k < npc-1; k+=2 ) { */
/* 		sum += T_n(i, xi)*sin(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	  }  */
/* 	} */
/*   } */

/*   /\* odd Chebyshev polynomials *\/ */
/*   for ( i = 1; i < nr-1; i += 2 ) { */
/* 	for ( j = 0; j < nt; j++ ) { */
/* 	  for ( k = 0; k < npc; k+=2 ) { */
/* 		sum += T_n(i, xi)*cos(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	  } */
/* 	  for ( k = 1; k < npc-1; k+=2 ) { */
/* 		sum += T_n(i, xi)*sin(j*theta)*(cos(k*phi) + sin(k*phi)); */
/* 	  }  */
/* 	} */
/*   } */
  
  return sum;
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

	/* transform for xi for fixed theta and phi */
	if(z == 0) { /* kernel */

	  /* even j ==> even Chebyshev series */
	  /* real part of coefficient */
	  for ( j = 0; j < nt; j += 2 ) {
		for ( k = 0; k < npc; k++ ) {
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
	  for ( j = 0; j < nt; j += 2 ) {
		/* k = 0 and k = npc-1 don't have immaginary parts */
		for ( k = 1; k < npc-1; k++ ) {
		  printf("in:   ");
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, IMAG );	
			printf("%f ", in1dxi[i]);
		  }
		  printf("\n");
		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  printf("out:  ");
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG, 
					   out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );	
			printf("%f ", coeff_get( coeff, z, i, j, k, IMAG ) );	
		  }  
		  printf("\n");
		}
	  } 
	  /* odd j ==> odd Chebyshev series */
	  /* real part of coefficient */
	  for ( j = 1; j < nt-1; j += 2 ) {
		for ( k = 1; k < npc; k++ ) {
		  for ( i = 1; i < nr; i++ ) {
			/* drop 1st point (it's zero) and reverse other points */
			in1dxiodd[i-1] = coeff_get( coeff, z, nr-i, j, k, REAL );
		  }
		  /* do type-3 DCT */
		  plan_forward = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT01, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  for ( i = 0; i < nr-1; i++ ) {
			coeff_set( coeff, z, i, j, k, REAL, 
					   out1dxiodd[i]/(nr-1) );
		  }
		  coeff_set( coeff, z, nr-1, j, k, REAL, 0.0 ); /* set coefficient for T_{nr}(xi) = 0 */
		}
	  }
	  /* odd j ==> odd Chebyshev series */
	  /* immaginary part of coefficient */
	  for ( j = 1; j < nt-1; j += 2 ) {
		/* k = 0 and k = npc-1 don't have immaginary parts */
		for ( k = 1; k < npc-1; k++ ) {
		  for ( i = 1; i < nr; i++ ) {
			/* drop 1st point (it's zero) and reverse other points */
			in1dxiodd[i-1] = coeff_get( coeff, z, nr-i, j, k, IMAG );
		  }
		  /* do type-3 DCT */
		  plan_forward = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT01, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  for ( i = 0; i < nr-1; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG, 
					   out1dxiodd[i]/(nr-1) );
		  }  
		  coeff_set( coeff, z, nr-1, j, k, IMAG, 0.0 ); /* set coefficient for T_{nr}(xi) = 0 */
		}
	  } 
	  
	} else { /* other zones */
	  
	  /* real part of coefficient */
	  for ( j = 0; j < nt; j++ ) {
		for ( k = 0; k < npc; k++ ) {
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, REAL );
		  }
		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, REAL, 
					   out1dcos[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
		  }
		}
	  }
	  /* immaginary part of coefficient */
	  for ( j = 0; j < nt; j++ ) {
		/* k = 0 and k = npc-1 don't have immaginary parts */
		for ( k = 1; k < npc-1; k++ ) {
		  for ( i = 0; i < nr; i++ ) {
			in1dxi[i] = coeff_get( coeff, z, i, j, k, IMAG );
		  }
		  plan_forward = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
		  fftw_execute ( plan_forward );
		  for ( i = 0; i < nr; i++ ) {
			coeff_set( coeff, z, i, j, k, IMAG, 
					   out1dcos[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
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
