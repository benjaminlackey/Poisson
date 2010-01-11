/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* fftw header */
#include <fftw3.h>

/* own header */
#include "poisson.h"


/*******************************************************************************************/
/* Finds surface where enthalpy drops to zero.                                             */
/* This function currently requires that the new surface be defined between zones 0 and 1. */
/*******************************************************************************************/
void findnewsurface(scalar2d *newsurface_scalar2d, scalar3d *enthalpy_scalar3d, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k;
  int nz, nr, nt, np;
  double enthalpy;
  double xi_i, xi_ip1, xi_ip2;
  double f_i, f_ip1, f_ip2;
  double xi, xiold;
  double f;
  double fprime;
  double acc = 1.0e-15;
  double alpha, beta;
  double rnew;

  nz = enthalpy_scalar3d->nz;
  nr = enthalpy_scalar3d->nr;
  nt = enthalpy_scalar3d->nt;
  np = enthalpy_scalar3d->np;
  
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {	
      
      /* find approximate position of surface */
      enthalpy = scalar3d_get(enthalpy_scalar3d, 1, 0, j, k); /* inside point of 1st shell */
      /*printf("enthalpy=%f\n", enthalpy);*/
      if(enthalpy >= 0.0) { /* H = 0 in zone 1 of old coordinate system */
	z = 1;
	i = 0;
	do {
	  enthalpy = scalar3d_get(enthalpy_scalar3d, z, i, j, k);
	  /*printf("enthalpy=%f\n", enthalpy);*/
	  i++;
	} while (enthalpy >= 0.0); /* i is just outside surface of star */
      } else {
	z = 0;
	i = nr-1; /* last radial point in kernel */
	do {
	  enthalpy = scalar3d_get(enthalpy_scalar3d, z, i, j, k);
	  i--;
	} while (enthalpy <= 0.0); /* i is just inside surface of star */
      }
      /*printf("z=%d\ti=%d\tj=%d\tk=%d\n", z, i, j, k);*/

      /* set (xi_i, f_i) points for interpolation */
      if (z == 0) {
	i = i-1;
	xi_i = sin(PI*i/(2*(nr-1)));
	xi_ip1 = sin(PI*(i+1)/(2*(nr-1)));
	xi_ip2 = sin(PI*(i+2)/(2*(nr-1)));
	f_i = scalar3d_get(enthalpy_scalar3d, z, i, j, k);
	f_ip1 = scalar3d_get(enthalpy_scalar3d, z, i+1, j, k);
	f_ip2 = scalar3d_get(enthalpy_scalar3d, z, i+2, j, k);      
	xi = xi_ip2; /* Initial guess for the zero */
      } else {
	i = i-2;
	xi_i = -cos(PI*i/(nr-1));
	xi_ip1 = -cos(PI*(i+1)/(nr-1));
	xi_ip2 = -cos(PI*(i+2)/(nr-1));
	f_i = scalar3d_get(enthalpy_scalar3d, z, i, j, k);
	f_ip1 = scalar3d_get(enthalpy_scalar3d, z, i+1, j, k);
	f_ip2 = scalar3d_get(enthalpy_scalar3d, z, i+2, j, k);
	xi = xi_i; /* Initial guess for the zero */
      }
      /*printf("z=%d\ti=%d\tj=%d\tk=%d\t%f\t%f\t%f\n", z, i, j, k, f_i, f_ip1, f_ip2);*/
      
      /* Start Newton's method to find the zero */
      do {
	/* Evaluate function at point xi.      */
	/* quadratic (3 point) interpolation */
	f = (xi - xi_ip1)*(xi - xi_ip2) / ((xi_i - xi_ip1)*(xi_i - xi_ip2))*f_i
	  + (xi - xi_i)*(xi - xi_ip2) / ((xi_ip1 - xi_i)*(xi_ip1 - xi_ip2))*f_ip1
	  + (xi - xi_i)*(xi - xi_ip1) / ((xi_ip2 - xi_i)*(xi_ip2 - xi_ip1))*f_ip2;	
	
	/* Evaluate derivative at point xi. */
	/* derivative of quadratic interpolation */
	fprime = ((xi - xi_ip1) + (xi - xi_ip2)) / ((xi_i - xi_ip1)*(xi_i - xi_ip2))*f_i
	  + ((xi - xi_i) + (xi - xi_ip2)) / ((xi_ip1 - xi_i)*(xi_ip1 - xi_ip2))*f_ip1
	  + ((xi - xi_i) + (xi - xi_ip1)) / ((xi_ip2 - xi_i)*(xi_ip2 - xi_ip1))*f_ip2;
	
	/* Make next guess for the zero */
	xiold = xi;
	xi = xiold - f / fprime;
      } while(fabs(xi-xiold) > acc); /* xi is in range [0, 1] or [-1, 1] */
      /*printf("z=%d\txi=%f\tj=%d\tk=%d\n", z, xi, j, k);*/      
      
      /* evaluate physical position r of surface */
      if(z == 0) {
	/* kernel */
	alpha = gsl_vector_get(alphalist, 0);
	rnew = alpha*(xi 
		      + xi*xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f_scalar2d, 0, j, k)
		      + 0.5*xi*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g_scalar2d, 0, j, k));
      } else { /* This breaks if star is in external zone */
	/* shells */
	alpha = gsl_vector_get(alphalist, z);
	beta = gsl_vector_get(betalist, z);
	rnew = alpha*(xi 
		      + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f_scalar2d, z, j, k)
		      + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g_scalar2d, z, j, k))
	  + beta;
      }
      /*printf("z=%d\tR=%f\tj=%d\tk=%d\n", z, rnew, j, k);*/
      scalar2d_set(newsurface_scalar2d, 0, j, k, rnew);
    }
  }
  
  /* just set new bound to old bound for z > 0 */
  /* this will probably change later */
  for ( z = 1; z < nz-1; z++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	alpha = gsl_vector_get(alphalist, z);
	beta = gsl_vector_get(betalist, z);
	rnew = alpha*(1.0 + scalar2d_get(g_scalar2d, z, j, k)) + beta;
	scalar2d_set(newsurface_scalar2d, z, j, k, rnew);
      }
    }
  }
  
}


/***************************************************/
/* This is an ideal place to parallelize the code. */
/* Each of nz*nr*nt*np evaluations is independent. */
/* This is also the most expensive operation.      */
/****************************************************************/
/* Take function evaluated at old gridpoints (funcold_scalar3d) */
/* and reevaluate them at new gridpoints (funcnew_scalar3d).    */
/****************************************************************/
void remapgrid(scalar3d *funcnew_scalar3d, gsl_vector *alphanew_vector, gsl_vector *betanew_vector, scalar2d *fnew_scalar2d, scalar2d *gnew_scalar2d, scalar3d *funcold_scalar3d, gsl_vector *alphaold_vector, gsl_vector *betaold_vector, scalar2d *fold_scalar2d, scalar2d *gold_scalar2d)
{
  int znew; /* zone of point in old coordinate system */
  int zold; /* zone of point in new coordinate system */
  int i, j, k;
  int nz, nr, nt, np;
  double xi, theta, phi;
  double rnew;
  double surfaceold;
  double alpha, beta;
  double f, g;
  double funcnew;
  scalar3d *rold_scalar3d, *rnew_scalar3d;
  coeff *funcold_coeff;

  nz = funcnew_scalar3d->nz;
  nr = funcnew_scalar3d->nr;
  nt = funcnew_scalar3d->nt;
  np = funcnew_scalar3d->np;
  
/*   printf("nz for fnew_scalar2d=%d\n", fnew_scalar2d->nz); */
/*   printf("nz for gnew_scalar2d=%d\n", gnew_scalar2d->nz); */
/*   printf("nz for fold_scalar2d=%d\n", fold_scalar2d->nz); */
/*   printf("nz for gold_scalar2d=%d\n", gold_scalar2d->nz); */

  rold_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  rnew_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  funcold_coeff = coeff_alloc(nz, nr, nt, np);

  /* find radius of old and new gridpoints */
  rofxtp(rold_scalar3d, alphaold_vector, betaold_vector, fold_scalar2d, gold_scalar2d);
  rofxtp(rnew_scalar3d, alphanew_vector, betanew_vector, fnew_scalar2d, gnew_scalar2d);
  
  /* evaluate coefficients of function in old coordinate system */
  gridtofourier(funcold_coeff, funcold_scalar3d, 0, 0);
  /*print_scalar3d(funcold_scalar3d);*/
  /*print_coeff(funcold_coeff);*/
  /* evaluate function at points (znew, rnew_i, theta_j, phi_k) in new coordinate system */
  for ( znew = 0; znew < nz-1; znew++ ) { /* deal with external zone later */
    /* MAY HAVE TO TREAT INNER OUTER POINTS (i=0, nr-1) SPECIALLY WHEN THERE ARE DISCONTINUITIES. */
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  rnew = scalar3d_get(rnew_scalar3d, znew, i, j, k); /* physical distance of new gridpoint */
	  /* find zold value of new gridpoint in old coordinate system */
	  zold = 0; /* start in kernel */
	  surfaceold = scalar3d_get(rold_scalar3d, zold, nr-1, j, k);
	  while ((rnew > surfaceold)&&(zold < nz-1)) {
	    zold++;
	    surfaceold = scalar3d_get(rold_scalar3d, zold, nr-1, j, k);
	  }
	  /*printf("(znew=%d, i=%d, j=%d, k=%d, zold=%d, nz=%d, nr=%d, nt=%d, np=%d)\n", znew, i, j, k, zold, nz, nr, nt, np);*/
	  if(zold==nz-1)
	    rnew = 1.0/rnew; /* now rnew = u = 1/r */
	  /*printf("(znew=%d, i=%d, j=%d, k=%d, zold=%d, nz=%d, nr=%d, nt=%d, np=%d)\n", znew, i, j, k, zold, nz, nr, nt, np);*/

	  /* find xi value of new gridpoint in old coordinate system */
	  alpha = gsl_vector_get(alphaold_vector, zold);
	  beta = gsl_vector_get(betaold_vector, zold);
	  f = scalar2d_get(fold_scalar2d, zold, j, k);
	  /*printf("hi\n");*/
	  g = scalar2d_get(gold_scalar2d, zold, j, k);
	  /*printf("hi2\n");*/
	  xi = xiofroru(nz, zold, alpha, beta, f, g, rnew);
	  /*printf("(znew=%d, i=%d, j=%d, k=%d, zold=%d, xiold=%.18e\n", znew, i, j, k, zold, xi);*/
	  theta =  PI*j/(nt-1);
	  phi = 2*PI*k/np;
	  /* do interpolation using coefficients in old coordinate system */
	  funcnew = spectralinterpolation(funcold_coeff, zold, xi, theta, phi);
	  /*printf("(znew=%d, i=%d, j=%d, k=%d, zold=%d, funcnew=%.18e\n", znew, i, j, k, zold, funcnew);*/
	  /* set function value at new gridpoint */
	  scalar3d_set(funcnew_scalar3d, znew, i, j, k, funcnew);
	  /*printf("(znew=%d, i=%d, j=%d, k=%d, zold=%d, R1=%f, rnew=%f, xiold=%f, theta=%f, phi=%f, f=%f)\n", znew, i, j, k, zold, scalar3d_get(rold_scalar3d, 1, 0, j, k), rnew, xi, theta, phi, funcnew);*/
	}
      }
    }
  }
  
  znew = nz-1; /* deal with external zone here */
  for ( i = 0; i < nr-1; i++ ) { /* deal with point at r = inf later */
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	rnew = 1/scalar3d_get(rnew_scalar3d, znew, i, j, k);
	/* find zold value of new gridpoint in old coordinate system */
	zold = 0; /* start in kernel */
	surfaceold = scalar3d_get(rold_scalar3d, zold, nr-1, j, k);
	while ((rnew > surfaceold)&&(zold < nz-1)) {
	  zold++;
	  surfaceold = scalar3d_get(rold_scalar3d, zold, nr-1, j, k);
	}
	if(zold==nz-1)
	  rnew = 1.0/rnew; /* now rnew = u = 1/r */
	
	/* find xi value of new gridpoint in old coordinate system */
	alpha = gsl_vector_get(alphaold_vector, zold);
	beta = gsl_vector_get(betaold_vector, zold);
	f = scalar2d_get(fold_scalar2d, zold, j, k);
	g = scalar2d_get(gold_scalar2d, zold, j, k);
	xi = xiofroru(nz, zold, alpha, beta, f, g, rnew);
	theta =  PI*j/(nt-1);
	phi = 2*PI*k/np;
	/* do interpolation using coefficients in old coordinate system */
	funcnew = spectralinterpolation(funcold_coeff, zold, xi, theta, phi);
	/* set function value at new gridpoint */
	scalar3d_set(funcnew_scalar3d, znew, i, j, k, funcnew);
      }
    }
  }
  
  /* deal with point at r = inf here */
  /* r = inf points in old, new coordinate system are the same */
  znew = nz-1;
  i = nr-1; 
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      scalar3d_set(funcnew_scalar3d, znew, i, j, k, scalar3d_get(funcold_scalar3d, znew, i, j, k));
    }
  }
  
  scalar3d_free(rold_scalar3d);
  scalar3d_free(rnew_scalar3d);
}

/******************************************************************************/
/* Root finder to find xi(r) given alpha, beta, f(theta, phi), g(theta, phi), */
/* where r(xi) is high order polynomial.                                      */
/* Currently using bisection which is slow (~50 function evaluations).        */
/******************************************************************************/
double xiofroru(int nz, int z, double alpha, double beta, double f, double g, double r)
{
  int i;
  int maxit = 100; /* maximum number of iterations */
  double xilow, ximid, xihigh;
  double ylow, ymid;
  double acc = 1.0e-15; /* absolute accuracy */

  /* set initial brackets */
  if (z == 0)
    xilow = 0.0;
  else
    xilow = -1.0;
  xihigh = 1.0;
  
  ylow = roruofxi(nz, z, alpha, beta, f, g, xilow) - r;
  for( i = 1; i <= maxit; i++ ) {
    ximid = 0.5*(xilow + xihigh); /* mid point */
    ymid = roruofxi(nz, z, alpha, beta, f, g, ximid) - r;
    
    if(ylow*ymid <= 0) /* root is bracketed by xilow and ximid */
      xihigh = ximid;
    else /* root is bracketed by ximid and xihigh */
      xilow = ximid;
    
    if(xihigh - xilow <= acc) {
      /*printf("%d\n", i);*/
      return ximid;
    }
  }

  /* too many iterations */
  printf("Number of iterations has exceeded %d in xiofroru.\n", maxit);
  printf("Accuracy is less than %e.\n", acc);
  return ximid;
}
  

/****************************************************************************/
/* Evaluate r(xi) or u(xi) given alpha, beta, f(theta, phi), g(theta, phi). */
/****************************************************************************/
double roruofxi(int nz, int z, double alpha, double beta, double f, double g, double xi)
{
  double r;

  if (z == 0) {
    r = alpha*(xi
	       + xi*xi*xi*xi*(3.0-2.0*xi*xi)*f
	       + 0.5*xi*xi*xi*(5.0-3.0*xi*xi)*g);
  } else if (z < nz-1) {
    r = alpha*(xi
	       + 0.25*(xi*xi*xi-3.0*xi+2.0)*f
	       + 0.25*(-xi*xi*xi+3.0*xi+2.0)*g)
      + beta;
  } else { /* external zone. r is actually u = 1/r here */
    r = alpha*(xi + 0.25*(xi*xi*xi-3.0*xi+2.0)*f - 1.0);
  }

  return r;
}


/*******************************************************************************/
/* Evaluate function, expressed as coefficients, at point (z, xi, theta, phi). */
/* Summing is done from highest order coefficients to smallest because         */
/* coefficients typically decrease exponentially with order.                   */
/*******************************************************************************/
double spectralinterpolation(coeff *func_coeff, int z, double xi, double theta, double phi)
{
  int i, j, k, imag;
  int nr, nt, np;
  double f, f_reim, f_k, f_jk;
  gsl_vector *coeff_vector, *coeffeven_vector, *coeffodd_vector;
  
  nr = func_coeff->nr;
  nt = func_coeff->nt;
  np = func_coeff->np;
  
  coeff_vector = gsl_vector_alloc(nr); /* {a_0, a_1, a_2, a_3, a_4, a_5, a_6} */
  coeffeven_vector = gsl_vector_alloc(2*nr-1); /* {a_0, 0, a_2, 0, a_4, 0, a_6} */
  coeffodd_vector = gsl_vector_alloc(2*nr-2); /* {0, a_1, 0, a_3, 0, a_5} */
  
  f = 0.0;
  for ( imag = 0; imag <= 1; imag++ ) {
    f_reim = 0.0;
    for ( k = np/2-imag; k >= imag; k-- ) { /* {0...np/2} for cos(k*phi).  {1, np/2-1} for sin(k*phi) */
      f_k = 0.0;
      for ( j = nt-1-k%2; j >= k%2; j-- ) { /* {0,...,nt-1} for k even. {1,...,nt-2} for k odd. */  
	/* sum over Chebyshev polynomials for each j, k */
	if (z > 0) { 
	  /* not in kernel */
	  for ( i = 0; i < nr; i++ ) {
	    gsl_vector_set(coeff_vector, i, coeff_get(func_coeff, z, i, j, k, imag));
	  }
	  f_jk = chebyshevinterpolation(coeff_vector, xi);
	} else if (!(j%2)) { 
	  /* in kernel and j is even */
	  gsl_vector_set_zero(coeffeven_vector);
	  for ( i = 0; i < nr; i++ ) {
	    gsl_vector_set(coeffeven_vector, 2*i, coeff_get(func_coeff, z, i, j, k, imag));
	  }
	  f_jk = chebyshevinterpolation(coeffeven_vector, xi);
	} else { 
	  /* in kernel and j is odd */
	  gsl_vector_set_zero(coeffodd_vector);
	  for ( i = 0; i < nr-1; i++ ) { /* last coefficient is even polynomial */
	    gsl_vector_set(coeffodd_vector, 2*i+1, coeff_get(func_coeff, z, i, j, k, imag));
	  }
	  /*printf("odd vector = ");
	    print_vector(coeffodd_vector);*/
	  f_jk = chebyshevinterpolation(coeffodd_vector, xi);
	  /*printf("xi=%.18e, f_jk=%.18e\n", xi, f_jk);*/
	}
	if (!(k%2)) { /* k is even */
	  f_k += f_jk*cos(j*theta);
	} else { /* k is odd */
	  f_k += f_jk*sin(j*theta);
	}
      }
      if (imag == 0) {
	f_reim += f_k*cos(k*phi);
      } else {
	f_reim += f_k*sin(k*phi);
      }
    }
    f += f_reim;
  }
  
  gsl_vector_free(coeff_vector);
  gsl_vector_free(coeffeven_vector);
  gsl_vector_free(coeffodd_vector);
  
  return f;
}

  
/****************************************************************/
/* Use Clenshaw's recurrence formula to evaluate function at xi */
/* from its Chebyshev coefficients.                             */
/****************************************************************/
double chebyshevinterpolation(gsl_vector *coeff_vector, double xi)
{ 
  int nr; /*number of polynomials*/
  int i;
  double y_i;
  double y_ip1; 
  double y_ip2;  
  
  nr = coeff_vector->size;
  
  /*Begin Clenshaw's recurrance formula.*/
  y_i = 0.0;
  y_ip1 = 0.0; 
  y_ip2 = 0.0;
  for (i = nr-1; i >= 1; i--) {
    y_ip2 = y_ip1;
    y_ip1 = y_i;
    y_i = 2*xi*y_ip1 - y_ip2 + gsl_vector_get(coeff_vector, i);
  }
  /*i is now 1*/
  return y_i*xi - y_ip1 + gsl_vector_get(coeff_vector, 0);
}
