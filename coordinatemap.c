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


/***********************************/
/* copy from source to destination */
/***********************************/
void scalar2d_memcpy(scalar2d *dest, scalar2d *source)
{
  int i;
  int nz, nt, np;

  nz = source->nz;
  nt = source->nt;
  np = source->np;

  for(i=0; i<nz*nt*np; i++)
    dest->data[i] = source->data[i];
}

void scalar3d_memcpy(scalar3d *dest, scalar3d *source)
{
  int i;
  int nz, nr, nt, np;

  nz = source->nz;  
  nr = source->nr;
  nt = source->nt;
  np = source->np;

  for(i=0; i<nz*nr*nt*np; i++)
    dest->data[i] = source->data[i];
}

void coeff_memcpy(coeff *dest, coeff *source)
{
  int i;
  int nz, nr, nt, np, npc;

  nz = source->nz;  
  nr = source->nr;
  nt = source->nt;
  np = source->np;

  npc = np/2 + 1;

  for(i=0; i<nz*nr*nt*2*npc; i++)
    dest->data[i] = source->data[i];
}

void bound_coeff_memcpy(bound_coeff *dest, bound_coeff *source)
{
  int i;
  int nz, nt, np, npc;

  nz = source->nz;  
  nt = source->nt;
  np = source->np;

  npc = np/2 + 1;

  for(i=0; i<nz*nt*2*npc; i++)
    dest->data[i] = source->data[i];
}
  
/*********************************/
/* Add a constant to a function. */
/*********************************/
void scalar3d_addconstant(scalar3d *in_grid, double number, scalar3d *sum_grid)
{
  int z, i, j, k;
  int nz, nr, nt, np;

  nz = sum_grid->nz;
  nr = sum_grid->nr;
  nt = sum_grid->nt;
  np = sum_grid->np;
  
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(sum_grid, z, i, j, k, 
		       scalar3d_get(in_grid, z, i, j, k) + number);
	}
      }
    }
  }  
}


/****************************************/
/* Adds 2 functions at the grid points. */
/****************************************/
void scalar3d_add(scalar3d *in1_grid, scalar3d *in2_grid, scalar3d *sum_grid)
{
  int z, i, j, k;
  int nz, nr, nt, np;

  nz = sum_grid->nz;
  nr = sum_grid->nr;
  nt = sum_grid->nt;
  np = sum_grid->np;
  
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(sum_grid, z, i, j, k, 
		       scalar3d_get(in1_grid, z, i, j, k) + scalar3d_get(in2_grid, z, i, j, k));
	}
      }
    }
  }  
}


/**********************************************/
/* Multiplies 2 functions at the grid points. */
/**********************************************/
void scalar3d_multiply(scalar3d *in1_grid, scalar3d *in2_grid, scalar3d *product_grid)
{
  int z, i, j, k;
  int nz, nr, nt, np;

  nz = product_grid->nz;
  nr = product_grid->nr;
  nt = product_grid->nt;
  np = product_grid->np;
  
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  scalar3d_set(product_grid, z, i, j, k, 
		       scalar3d_get(in1_grid, z, i, j, k) * scalar3d_get(in2_grid, z, i, j, k));
	}
      }
    }
  }  
}


/*******************************************************************************************************/
/*  Combine the functions map_physicaltogrid_kernel, map_physicaltogrid_shell, map_physicaltogrid_ext  */
/*******************************************************************************************************/
void map_physicaltogrid(scalar2d *b_zjk, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f_grid, scalar2d *g_grid)
{
  int nz;
  int z;
  
  nz = b_zjk->nz + 1;

  map_physicaltogrid_kernel(b_zjk, alphalist, f_grid, g_grid);
  for(z=1; z<nz-1; z++) {
    map_physicaltogrid_shell(b_zjk, z, alphalist, betalist, f_grid, g_grid);
  }
  map_physicaltogrid_ext(b_zjk, alphalist, f_grid);
  
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
  
/*   /\* Initialize these structures to 0 *\/ */
/*   /\* you should probably just make an initialization function instead *\/ */
/*   for(z=0; z<nz-1; z++) { */
/*     for(j=0; j<nt; j++) { */
/*       for(k=0; k<npc; k++) { */
/* 	/\*printf("z=%d, j=%d, k=%d\n", z, j, k);*\/ */
/* 	scalar2d_set(f_grid, z, j, k, 0.0); */
/* 	bound_coeff_set(f_fourier, z, j, k, REAL, 0.0); */
/* 	bound_coeff_set(f_fourier, z, j, k, IMAG, 0.0); */
/*       } */
/*     } */
/*   } */
/*   for(z=0; z<nz; z++) { */
/*     for(j=0; j<nt; j++) { */
/*       for(k=0; k<npc; k++) { */
/* 	/\*printf("z=%d, j=%d, k=%d\n", z, j, k);*\/ */
/* 	bound_coeff_set(fodd_fourier, z, j, k, REAL, 0.0); */
/* 	bound_coeff_set(fodd_fourier, z, j, k, IMAG, 0.0); */
/* 	bound_coeff_set(geven_fourier, z, j, k, REAL, 0.0); */
/* 	bound_coeff_set(geven_fourier, z, j, k, IMAG, 0.0); */
/*       } */
/*     } */
/*   } */
  
  /* evaluate some special point and find f = S - S_point for each point on kernel boundary */
  specialpoint = scalar2d_get(b_zjk, 0, 0, 0);
  /*printf("specialpoint=%f\n", specialpoint);*/
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
	  scalar2d_set( f_grid, 0, j, k, 
			scalar2d_get(b_zjk, 0, j, k) - specialpoint );
    }
  }
  /*printf("f_grid = b_zjk - specialpoint\n");*/
  /*print_scalar2d(f_grid);*/

  /* decompose f into double fourier series */
  gridtofourier_bound(f_fourier, f_grid);
  /*print_bound_coeff(f_fourier);*/

  /* set fodd~ to odd phi harmonics */
  for(k=1; k<npc; k += 2) {
    for(j=0; j<nt; j++) {
      bound_coeff_set(fodd_fourier, 0, j, k, REAL, bound_coeff_get(f_fourier, 0, j, k, REAL));
      bound_coeff_set(fodd_fourier, 0, j, k, IMAG, bound_coeff_get(f_fourier, 0, j, k, IMAG));
    }
  }
  /*printf("fodd_fourier:\n");
    print_bound_coeff(fodd_fourier);*/

  /* set geven~ to even phi harmonics */
  for(k=0; k<npc; k += 2) {
    for(j=0; j<nt; j++) {
      bound_coeff_set(geven_fourier, 0, j, k, REAL, bound_coeff_get(f_fourier, 0, j, k, REAL));
      bound_coeff_set(geven_fourier, 0, j, k, IMAG, bound_coeff_get(f_fourier, 0, j, k, IMAG));
    }
  }
  /*printf("geven_fourier:\n");
    print_bound_coeff(geven_fourier);*/  

  
  /* transform fodd~ and geven~ back to spatial domain */
  fouriertogrid_bound(fodd_grid, fodd_fourier, 0);
  fouriertogrid_bound(geven_grid, geven_fourier, 0);
  /*printf("fodd_grid:\n");
  print_scalar2d(fodd_grid);
  printf("fodd_grid:\n");
  print_scalar2d(fodd_grid);*/
  
  /* calculate mu (the negative of the minimum value of geven~) */
  mu = scalar2d_get(geven_grid, 0, 0, 0);
  for(k=0; k<np; k++) {
    for(j=0; j<nt; j++) {
      new = scalar2d_get(geven_grid, 0, j, k);
      mu = (new < mu ? new : mu);
    }
  }
  mu = -mu;
  /*printf("mu=%f\n", mu);*/

  /* calculate alpha */
  alpha = specialpoint - mu;
  gsl_vector_set(alphalist, 0, alpha);
  /*printf("alpha=%f\n", alpha);*/

  /* rescale fodd~ and geven~ to calculate fodd and geven */ 
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
      /* fodd: */
      scalar2d_set(fodd_grid, 0, j, k, scalar2d_get(fodd_grid, 0, j, k)/alpha);
      /* geven: */
      scalar2d_set(geven_grid, 0, j, k, (scalar2d_get(geven_grid, 0, j, k) + mu) / alpha);
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
  specialpointout = scalar2d_get(b_zjk, z, 0, 0);

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
  /*printf("specialpointin=%f\n", specialpointin);*/
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
  /*printf("lambda=%f\n", lambda);*/

  /* calculate alpha */
  alpha = 0.5*(lambda - 1.0/specialpointin);
  gsl_vector_set(alphalist, nz-1, alpha);
  /*printf("alpha=%f\n", alpha);*/

  /* rescale f~ to calculate f */ 
  for(k=0; k<np; k++) { 
    for(j=0; j<nt; j++) {
      scalar2d_set(f_grid, nz-1, j, k, (scalar2d_get(f_grid, nz-1, j, k) + lambda)/alpha);
    }
  }
  
}


/*********************************************************/
/* Should be called ratgridpoints or something like that */
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


/******************************************************************************************************/
/* Take a function for the boundary b = b(z, theta, phi) and evaluate it on the boundary grid points. */
/******************************************************************************************************/ 
void boundarytogrid(scalar2d *boundary_scalar2d, double (*boundary)(int z, double theta, double phi))
{
  int z, j, k;
  int nz, nt, np;
  double theta_j, phi_k;
  
  nz = boundary_scalar2d->nz + 1;
  nt = boundary_scalar2d->nt;
  np = boundary_scalar2d->np;
  
  for ( z = 0; z < nz-1; z++ ) {
    for ( k = 0; k < np; k++ ) {
      phi_k = 2*PI*k/np;
      for ( j = 0; j < nt; j++ ) {
	theta_j = PI*j/(nt-1);
	scalar2d_set(boundary_scalar2d, z, j, k, boundary(z, theta_j, phi_k));
      }
    }
  }
}


/******************************************************************************************************/
/* Take a function f = f(z, r, theta, phi) and evaluate it on the surface matched grid func_scalar3d. */
/******************************************************************************************************/ 
void functiontogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi))
{
  int z, i, j, k;
  int nz, nr, nt, np;
  double xi_i, theta_j, phi_k;
  double r;
  scalar3d *r_scalar3d;

  nz = func_scalar3d->nz;
  nr = func_scalar3d->nr;
  nt = func_scalar3d->nt;
  np = func_scalar3d->np;
  
  /* value of R (or U) at each point (xi, theta, phi) on grid */
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);  
  
  /* Find the radial position of each point. */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);

  /* kernel and shells */
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
	  r = scalar3d_get(r_scalar3d, z, i, j, k);
	  scalar3d_set(func_scalar3d, z, i, j, k, func(z, r, theta_j, phi_k));
	}
      }
    }
  }
  /* external domain except for r = infinity */
  z=nz-1;
  for ( i = 0; i < nr-1; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	xi_i = -cos(PI*i/(nr-1));
	theta_j = PI*j/(nt-1);
	phi_k = 2*PI*k/np;
	r = scalar3d_get(r_scalar3d, z, i, j, k);
	scalar3d_set(func_scalar3d, z, i, j, k, func(z, 1/r, theta_j, phi_k));
      }
    }
  }
  /* set to zero at r = infinity */
  z=nz-1;
  i=nr-1;
  for ( j = 0; j < nt; j++ ) {
    for ( k = 0; k < np; k++ ) {
      scalar3d_set(func_scalar3d, z, i, j, k, 0.0);
    }
  }

  scalar3d_free(r_scalar3d);
}


/***************************************************************************************/
/* Take a function f = f(z, xi, theta, phi) and evaluate it on the grid func_scalar3d. */
/***************************************************************************************/ 
void functiontogrid_xi(scalar3d *func_scalar3d, double (*func)(int z, double xi, double theta, double phi))
{
  int z, i, j, k;
  int nz, nr, nt, np;
  double xi_i, theta_j, phi_k;
  
  nz = func_scalar3d->nz;
  nr = func_scalar3d->nr;
  nt = func_scalar3d->nt;
  np = func_scalar3d->np;
  
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  xi_i = ((z==0) ? sin(PI*i/(2.0*(nr-1))) : -cos(PI*i/(nr-1)));
	  theta_j = PI*j/(nt-1);
	  phi_k = 2*PI*k/np;
	  scalar3d_set(func_scalar3d, z, i, j, k, func(z, xi_i, theta_j, phi_k));
	}
      }
    }
  }
}

/* /\*************************************************************\/ */
/* /\* Use recursion relations to divide f in kernel by xi.      *\/ */
/* /\* This is used to analytically get rid of the zero xi=R=0.  *\/ */
/* /\*************************************************************\/ */
/* dividebyxi_kernel(coeff fbyxi_coeff, coeff f_coeff) */
/* { */
/*   int z; */
/*   int i; */
/*   int j; */
/*   int k; */
/*   int nz; */
/*   int nr; */
/*   int nt; */
/*   int np; */
/*   double fbyxi_i; */
/*   double f_ip1; */
/*   double fbyxi_ip2; */
  
/*   nz = f_coeff->nz; */
/*   nr = f_coeff->nr; */
/*   nt = f_coeff->nt; */
/*   np = f_coeff->np; */

/*   /\* even kernel *\/ */
/*   /\* even Chebyshev polynomials become odd Chebyshev polynomials *\/ */
/*   /\* even Chebyshev polynomials are indexed as T_{2i}(xi) *\/ */
/*   /\* odd Chebyshev polynomials are indexed as T_{2i+1}(xi) *\/ */
/*   z = 0; */
/*   for(imag=0; imag<=1; imag++) { */
/*     for(k = imag; k < npc-imag; k++) { /\* imaginary parts for k=0 and k=npc-1 are zero *\/ */
/*       for(j = 0; j < nt; j+=2) { /\* even j only *\/ */
/* 	/\* do recursion for 1/xi *\/ */
/* 	/\* set 1/xi coefficient for T_{2(nr-1)+1} to zero: *\/ */
/* 	coeff_set(fbyxi, z, nr-1, j, k, imag, 0.0); */
/* 	/\* set fbyxi coefficient for T_{2(nr-2)+1}: *\/ */
/* 	fbyxi_i = 2.0*coeff_get(f, z, nr-1, j, k, imag); */
/* 	coeff_set(fbyxi, z, nr-2, j, k, imag, fprime_i); */
/* 	/\* now set fbyxi for nr-3 and below: *\/ */
/* 	fbyxi_ip2 = 0.0; /\* where i starts at nr-3 below *\/ */
/* 	for(i=nr-3; i>=0; i--) { */
/* 	  f_ip1 = coeff_get(f, z, i+1, j, k, imag); */
/* 	  fbyxi_ip2 = coeff_get(fbyxi, z, i+1, j, k, imag); */
/* 	  fbyxi_i = 2.0*f_ip1 + fbyxi_ip2; */
/* 	  coeff_set(fbyxi, z, i, j, k, imag, fbyxi_i); */
/* 	}  */
/*       } */
/*     } */
/*   } */

/*   /\* odd kernel *\/ */
/*   /\* odd Chebyshev polynomials become even Chebyshev polynomials *\/ */
/*   for(imag=0; imag<=1; imag++) { */
/*     for(k = imag; k < npc-imag; k++) { /\* imaginary parts for k=0 and k=npc-1 are zero *\/ */
/*       for(j = 1; j < nt-1; j+=2) { /\* odd j only *\/ */
/* 	/\* do recursion for 1/xi *\/ */
/* 	/\* set 1/xi coefficient for T_{2(nr-1)} to zero: *\/ */
/* 	coeff_set(fbyxi, z, nr-1, j, k, imag, 0.0); */
/* 	/\* set fbyxi coefficient for T_{2(nr-2)}: *\/ */
/* 	fbyxi_i = 2.0*coeff_get(f, z, nr-2, j, k, imag); */
/* 	coeff_set(fbyxi, z, nr-2, j, k, imag, fprime_i); */
/* 	/\* now set dfdxi for nr-3 and below: *\/ */
/* 	fbyxi_ip2 = 0.0; /\* where i starts at nr-3 below *\/ */
/* 	for(i=nr-3; i>=0; i--) { */
/* 	  f_ip1 = coeff_get(f, z, i, j, k, imag); */
/* 	  fbyxi_ip2 = coeff_get(fbyxi, z, i+1, j, k, imag); */
/* 	  fbyxi_i = 2*f_ip1 + fbyxi_ip2; */
/* 	  coeff_set(fbyxi, z, i, j, k, imag, fprime_i); */
/* 	}  */
/*       } */
/*     } */
/*   } */

/* } */


/* /\**********************************************************************\/ */
/* /\* Use recursion relations to divide f in external domain by (xi-1).  *\/ */
/* /\* This is used to analytically get rid of the zero xi=1, U=0, R=inf. *\/ */
/* /\* This function is not done!!!!!!!!!!!!!!                            *\/ */
/* /\**********************************************************************\/ */
/* dividebyximin1_ext(coeff fbyxim1_coeff, coeff f_coeff) */
/* { */
/*   int z; */
/*   int i; */
/*   int j; */
/*   int k; */
/*   int nz; */
/*   int nr; */
/*   int nt; */
/*   int np; */
/*   double fbyxi_i; */
/*   double f_ip1; */
/*   double fbyxi_ip2; */
  
/*   nz = f_coeff->nz; */
/*   nr = f_coeff->nr; */
/*   nt = f_coeff->nt; */
/*   np = f_coeff->np; */

/*   z = nz-1; */
/*   for(imag=0; imag<=1; imag++) { */
/*     for(k = imag; k < npc-imag; k++) { /\* imaginary parts for k=0 and k=npc-1 are zero *\/ */
/*       for(j = 0; j < nt; j++) { */
/* 	/\* do recursion for 1/(xi-1) *\/ */
/* 	/\* set 1/(xi-1) coefficient for T_{nr-1} to zero: *\/ */
/* 	coeff_set(fbyxim1_coeff, z, nr-1, i, j, k, imag, 0.0); */
/* 	/\* set fbyxi coefficient for T_{2(nr-2)+1}: *\/ */
/* 	fbyxi_i = 2.0*coeff_get(f, z, nr-1, j, k, imag); */
/* 	coeff_set(fbyxi, z, nr-2, j, k, imag, fprime_i); */
/* 	/\* now set fbyxi for nr-3 and below: *\/ */
/* 	fbyxi_ip2 = 0.0; /\* where i starts at nr-3 below *\/ */
/* 	for(i=nr-3; i>=0; i--) { */
/* 	  f_ip1 = coeff_get(f, z, i+1, j, k, imag); */
/* 	  fbyxi_ip2 = coeff_get(fbyxi, z, i+1, j, k, imag); */
/* 	  fbyxi_i = 2.0*f_ip1 + fbyxi_ip2; */
/* 	  coeff_set(fbyxi, z, i, j, k, imag, fbyxi_i); */
/* 	}  */
/*       } */
/*     } */
/*   } */

/* } */

/********************************************************************************************/
/* Take function in coefficient space and divide by R(xi, theta, phi) or U(xi, theta, phi). */
/* Returns the value at the surface-matched collocation points.                             */
/* r=0 and r=inf need to be treated analytically.                                           */
/********************************************************************************************/
void dividebyr(scalar3d *fbyr_scalar3d, coeff *f_coeff, int xishift, int thetashift, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g)
{
  int z, i, j, k, imag;
  int nz, nr, nt, np, npc;
  double alpha;
  double beta;
  double xi;
  double f_gridpoint;
  double den;
  double sum;
  bound_coeff *fbyxi_bound_coeff;
  scalar2d *fbyxi_scalar2d;

  nz = f_coeff->nz;
  nr = f_coeff->nr;
  nt = f_coeff->nt;
  np = f_coeff->np;

  npc = ( np / 2 ) + 1;  

  fbyxi_bound_coeff = bound_coeff_alloc(1, nt, np);
  fbyxi_scalar2d = scalar2d_alloc(1, nt, np);

  /*>>>>>>>>>>> EVALUATE F/R AT R=0 <<<<<<<<<<<<<<*/
  
  /* Divide by xi at xi = 0. */
  /* Use l'Hopital's rule to take lim_{xi->0} f(xi)/xi. */
  /* Only odd Chebyshev polynomials survive. */
  /*     -odd Chebyshev series if j is odd or */
  /*     -j is even and coeff was acted on by an odd number of d/dxi derivitives */  
  for(imag=0; imag<=1; imag++) {
    for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */ 
      for(j = 1-xishift; j < nt-1+xishift; j += 2) {
	sum = 0.0;
	for(i=nr-2+xishift; i>=0; i--) {
	  sum += neg1toi(i)*(2*i+1)*coeff_get(f_coeff, 0, i, j, k, imag);
	}
	bound_coeff_set(fbyxi_bound_coeff, 0, j, k, imag, sum);
      }
    }
  }
  print_bound_coeff(fbyxi_bound_coeff);

  fouriertogrid_bound(fbyxi_scalar2d, fbyxi_bound_coeff, thetashift);
  print_scalar2d(fbyxi_scalar2d); 

  /* convert to values on grid: */
  /* fbyr_grid is still just f at grid points */
  fouriertogrid(fbyr_scalar3d, f_coeff, xishift, thetashift);

  /* evaluate at r=0 for each angle theta, phi */
  alpha = gsl_vector_get(alphalist, 0);
  for(j=0; j<nt; j++) {
    for(k=0; k<np; k++) {
      scalar3d_set(fbyr_scalar3d, 0, 0, j, k, scalar2d_get(fbyxi_scalar2d, 0, j, k) / alpha);
      /*printf("k=%d, j=%d, sum=%f\n", k, j, scalar3d_get(fbyr_scalar3d, 0, 0, j, k));*/
    }
  }

  /* print_coeff(f_coeff); */
/*   /\* Evaluate f/R at r = 0. *\/ */
/*   /\* set theta = phi = 0 and use l'Hopital's rule to take lim_{xi->0} f'(xi)/R'(xi) *\/ */
/*   /\* the limit should be the same from any direction (theta, phi). *\/ */
/*   /\* start adding from highest order coefficients (they have smallest values) *\/ */
/*   alpha = gsl_vector_get(alphalist, 0); */
/*   sum = 0.0; */
/*   for(k=thetashift; k<npc; k+=2) { */
/*     for(j=1-xishift; j<nt-1+xishift; j+=2) { */
/*       for(i=0; i<nr-1+xishift; i++) { */
/* 	sum += neg1toi(i)*(2*i+1)*coeff_get(f_coeff, 0, i, j, k, REAL); */
/* 	printf("k=%d, j=%d, i=%d, coeff=%.18e term=%.18e\n", k, j, i, coeff_get(f_coeff, 0, i, j, k, REAL), neg1toi(i)*(2*i+1)*coeff_get(f_coeff, 0, i, j, k, REAL)); */
/*       } */
/*     } */
/*   } */
/*   printf("sum = %f\n", sum); */
/*   printf("alpha = %f\n", alpha);    */
/*   sum /= alpha; */
  


/*   /\* convert to values on grid: *\/ */
/*   /\* fbyr_grid is still just f at grid points *\/ */
/*   fouriertogrid(fbyr_grid, f_coeff, xishift, thetashift); */

/*   /\* should be same for all angles *\/ */
/*   for(j=0; j<nt; j++) { */
/*     for(k=0; k<np; k++) { */
/*       scalar3d_set(fbyr_grid, 0, 0, j, k, sum); */
/*       /\*printf("k=%d, j=%d, sum=%f\n", k, j, scalar3d_get(fbyr_grid, 0, 0, j, k));*\/ */
/*     } */
/*   }   */

 /*>>>>>>>>>>> EVALUATE F/R EVERYWHERE ELSE <<<<<<<<<<<<<<*/

  /* kernel (except for r = 0) */
  /* return result in same structure */
  for(i=1; i<nr; i++) { /* Treat the i = 0 (r = 0) case separately (above). */
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	f_gridpoint = scalar3d_get(fbyr_scalar3d, 0, i, j, k);
	den = alpha*(xi
		     + xi*xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f, 0, j, k)
		     + 0.5*xi*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g, 0, j, k));
	scalar3d_set(fbyr_scalar3d, 0, i, j, k, f_gridpoint/den);
      }
    }
  }

  /* shells */
  /* R has no zeros in the shells */
  for(z=1; z<nz-1; z++) {
    
    alpha = gsl_vector_get(alphalist, z);
    beta = gsl_vector_get(betalist, z);
    for(i=0; i<nr; i++) {
      xi = -cos(PI*i/(nr-1));
      for(j=0; j<nt; j++) {
	for(k=0; k<np; k++) {
	  f_gridpoint = scalar3d_get(fbyr_scalar3d, z, i, j, k);
	  den = alpha*(xi
		       + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, z, j, k)
		       + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g, z, j, k))
	    + beta;
	  scalar3d_set(fbyr_scalar3d, z, i, j, k, f_gridpoint/den);
	}
      }
    }
    
  }
  
  /* external domain */
  alpha = gsl_vector_get(alphalist, nz-1);
  for(i=0; i<nr-1; i++) { /* Treat the i = nr-1 (r = inf) case separately. */
    xi = -cos(PI*i/(nr-1));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	f_gridpoint = scalar3d_get(fbyr_scalar3d, nz-1, i, j, k);
	den = alpha*(xi + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, nz-1, j, k) - 1.0);
	scalar3d_set(fbyr_scalar3d, nz-1, i, j, k, f_gridpoint/den);
      }
    }
  }
  /* set f/U to 0.0 at r = inf since f dacays faster than 1/r */
  for(j=0; j<nt; j++) {
    for(k=0; k<np; k++) {
      scalar3d_set(fbyr_scalar3d, nz-1, nr-1, j, k, 0.0);
    }
  }
  
  bound_coeff_free(fbyxi_bound_coeff);
  scalar2d_free(fbyxi_scalar2d);
}


/**************************************************************/
/* J_1 = dR(xi, theta, phi)/dxi                               */
/* Evaluate J_1 at surface-matched coordinate points.         */
/**************************************************************/
void jacobian1(scalar3d *jacobian, gsl_vector *alphalist, scalar2d *f, scalar2d *g)
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
  double xi;
  double jpoint; /* jacobian J_1 evaluated at a specific point */
  
  nz = jacobian->nz;
  nr = jacobian->nr;
  nt = jacobian->nt;
  np = jacobian->np;
  
  /* kernel */
  alpha = gsl_vector_get(alphalist, 0);
  for(i=0; i<nr; i++) {
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	jpoint = alpha*(1.0
			+ 12.0*xi*xi*xi*(1.0-xi*xi)*scalar2d_get(f, 0, j, k)
			+ 7.5*xi*xi*(1.0-xi*xi)*scalar2d_get(g, 0, j, k));
	scalar3d_set(jacobian, 0, i, j, k, jpoint);
      }
    }
  }
  
  /* shells */
  for(z=1; z<nz-1; z++) {
    
    alpha = gsl_vector_get(alphalist, z);
    for(i=0; i<nr; i++) {
      xi = -cos(PI*i/(nr-1));
      for(j=0; j<nt; j++) {
	for(k=0; k<np; k++) {
	  jpoint = alpha*(1.0
			  + 0.75*(xi*xi-1.0)*scalar2d_get(f, z, j, k)
			  + 0.75*(1.0-xi*xi)*scalar2d_get(g, z, j, k));
	  scalar3d_set(jacobian, z, i, j, k, jpoint);
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
	jpoint = alpha*(1.0 + 0.75*(xi*xi-1.0)*scalar2d_get(f, nz-1, j, k));
	scalar3d_set(jacobian, nz-1, i, j, k, jpoint);
      }
    }
  }
  
}


/**********************************/
/* d^2 R(xi, theta, phi) / d xi^2 */
/**********************************/
void d2rdxi2(scalar3d *rxx, gsl_vector *alphalist, scalar2d *f, scalar2d *g)
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
  double xi;
  double rxxpoint; /* value of R_{xi xi} evaluated at a specific point */
  
  nz = rxx->nz;
  nr = rxx->nr;
  nt = rxx->nt;
  np = rxx->np;
  
  /* kernel */
  alpha = gsl_vector_get(alphalist, 0);
  for(i=0; i<nr; i++) {
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	rxxpoint = alpha*(xi*xi*(36.0-60.0*xi*xi)*scalar2d_get(f, 0, j, k)
			  + xi*(15.0-30.0*xi*xi)*scalar2d_get(g, 0, j, k));
	scalar3d_set(rxx, 0, i, j, k, rxxpoint);
      }
    }
  }
  
  /* shells */
  for(z=1; z<nz-1; z++) {
    
    alpha = gsl_vector_get(alphalist, z);
    for(i=0; i<nr; i++) {
      xi = -cos(PI*i/(nr-1));
      for(j=0; j<nt; j++) {
	for(k=0; k<np; k++) {
	  rxxpoint = 1.5*alpha*xi*(scalar2d_get(f, z, j, k) - scalar2d_get(g, z, j, k));
	  scalar3d_set(rxx, z, i, j, k, rxxpoint);
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
	rxxpoint = 1.5*alpha*xi*scalar2d_get(f, nz-1, j, k);
	scalar3d_set(rxx, nz-1, i, j, k, rxxpoint);
      }
    }
  }
  
}

/***********************************************/
/* Take partial derivative with respect to xi. */
/*           d f / d xi                        */
/***********************************************/
void dfdxi(coeff *dfdxi, coeff *f)
{
  int z;
  int i;
  int j;
  int k;
  int imag; /* imag=0 is real.  imag=1 is imaginary. */
  int nz;
  int nr; 
  int nt;
  int np;
  int npc;
  double fprime_i; /* f'_i(xi) */
  double f_ip1; /* f_{i+1}(xi) */
  double fprime_ip2; /* f'_{i+3}(xi) */
  
  nz = f->nz;
  nr = f->nr;
  nt = f->nt;
  np = f->np;
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* even kernel */
  /* even Chebyshev polynomials become odd Chebyshev polynomials */
  /* even Chebyshev polynomials are indexed as T_{2i}(xi) */
  /* odd Chebyshev polynomials are indexed as T_{2i+1}(xi) */
  z = 0;
  for(imag=0; imag<=1; imag++) {
    for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */
      for(j = 0; j < nt; j+=2) { /* even j only */
	/* do recursion for xi derivatives */
	/* set dfdxi coefficient for T_{2(nr-1)+1} to zero: */
	coeff_set(dfdxi, z, nr-1, j, k, imag, 0.0);
	/* set dfdxi coefficient for T_{2(nr-2)}: */
	fprime_i = 4*(nr-1)*coeff_get(f, z, nr-1, j, k, imag);
	coeff_set(dfdxi, z, nr-2, j, k, imag, fprime_i);
	/* now set dfdxi for nr-3 and below: */
	fprime_ip2 = 0.0; /* where i starts at nr-3 below */
	for(i=nr-3; i>=0; i--) {
	  f_ip1 = coeff_get(f, z, i+1, j, k, imag);
	  fprime_ip2 = coeff_get(dfdxi, z, i+1, j, k, imag);
	  fprime_i = 4*(i+1)*f_ip1 + fprime_ip2;
	  coeff_set(dfdxi, z, i, j, k, imag, fprime_i);
	} 
      }
    }
  }
  
  /* odd kernel */
  /* odd Chebyshev polynomials become even Chebyshev polynomials */
  for(imag=0; imag<=1; imag++) {
    for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */
      for(j = 1; j < nt-1; j+=2) { /* odd j only */
	/* do recursion for xi derivatives */
	/* set dfdxi coefficient for T_{2(nr-1)} to zero: */
	coeff_set(dfdxi, z, nr-1, j, k, imag, 0.0);
	/* set dfdxi coefficient for T_{2(nr-2)+1}: */
	fprime_i = 2*(2*nr-3)*coeff_get(f, z, nr-2, j, k, imag);
	coeff_set(dfdxi, z, nr-2, j, k, imag, fprime_i);
	/* now set dfdxi for nr-3 and below: */
	fprime_ip2 = 0.0; /* where i starts at nr-3 below */
	for(i=nr-3; i>=0; i--) {
	  f_ip1 = coeff_get(f, z, i, j, k, imag);
	  fprime_ip2 = coeff_get(dfdxi, z, i+1, j, k, imag);
	  fprime_i = (2*(2*i+1)*f_ip1 + fprime_ip2)/(1+delta(i, 0));
	  coeff_set(dfdxi, z, i, j, k, imag, fprime_i);
	} 
      }
    }
  }
  
  
  /* shells and external domain */
  for(z=1; z<nz; z++) {
    for(imag=0; imag<=1; imag++) {
      for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j = 0; j < nt; j++) { /* all j */
	  /* do recursion for xi derivatives */
	  /* set dfdxi for nr-1: */
	  coeff_set(dfdxi, z, nr-1, j, k, imag, 0.0);
	  /* set dfdxi for nr-2: */
	  fprime_i = 2*(nr-1)*coeff_get(f, z, nr-1, j, k, imag);
	  coeff_set(dfdxi, z, nr-2, j, k, imag, fprime_i);
	  /* now set dfdxi for nr-3 and below: */
	  fprime_ip2 = 0.0; /* where i starts at nr-3 below, so f'_{i+2} = f'_{nr-1} = 0 */
	  for(i=nr-3; i>=0; i--) {
	    f_ip1 = coeff_get(f, z, i+1, j, k, imag);
	    fprime_ip2 = coeff_get(dfdxi, z, i+2, j, k, imag);
	    fprime_i = (2*(i+1)*f_ip1 + fprime_ip2)/(1+delta(i, 0));
	    coeff_set(dfdxi, z, i, j, k, imag, fprime_i);
	  } 
	}
      }
    }
  }
  
  /* explicitly set the components of dfdxi that are 0 to 0 just in case dfdxi wasn't initialized ? */
}


/**************************************************************/
/* J_2 = 1/R * dR(xi, theta, phi)/dtheta                      */
/**************************************************************/
void jacobian2(scalar3d *jacobian, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g)
{
  int imag;
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  double alpha;
  double beta;
  double xi;
  double num;
  double den;
  bound_coeff *f_fourier;
  bound_coeff *g_fourier;
  scalar2d *dfdtheta;
  scalar2d *dgdtheta;
  
  nz = jacobian->nz;
  nr = jacobian->nr;
  nt = jacobian->nt;
  np = jacobian->np;
  
  npc = ( np / 2 ) + 1;

  /* allocate memory here */
  f_fourier = bound_coeff_alloc(nz, nt, np);
  g_fourier = bound_coeff_alloc(nz, nt, np);
  dfdtheta = scalar2d_alloc(nz, nt, np);
  dgdtheta = scalar2d_alloc(nz, nt, np);
  
  /* go to Fourier space */
  gridtofourier_bound(f_fourier, f);
  gridtofourier_bound(g_fourier, g);
  
  /* take d/dtheta derivatives of f and g for all zones */
  /* and return them in the same structure */
  
  /* d/dt[a_j cos(j theta)] = -j a_j sin(j theta) */  
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) {   
      for(k=2*imag; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j=0; j<nt; j++) {
	  /* f: */
	  bound_coeff_set(f_fourier, z, j, k, imag, 
			 -j * bound_coeff_get(f_fourier, z, j, k, imag));
	  if(z < nz-1) { /* g is not defined in external zone (nz-1) */
	    /* g: */
	    bound_coeff_set(g_fourier, z, j, k, imag, 
			   -j * bound_coeff_get(g_fourier, z, j, k, imag));
	  }
	}
      }
    }
  }

  /* d/dt[a_j sin(j theta)] = +j a_j cos(j theta) */
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) {   
      for(k=1; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j=1; j<nt-1; j++) {  
	  /* f: */
	  bound_coeff_set(f_fourier, z, j, k, imag, 
			 j * bound_coeff_get(f_fourier, z, j, k, imag));
	  if(z < nz-1) { /* g is not defined in external zone (nz-1) */
	    /* g: */
	    bound_coeff_set(g_fourier, z, j, k, imag, 
			   j * bound_coeff_get(g_fourier, z, j, k, imag));
	  }
	}
      }
    }
  }
  
  /* return from Fourier space to grid */
  fouriertogrid_bound(dfdtheta, f_fourier, 1);
  fouriertogrid_bound(dgdtheta, g_fourier, 1);
  
  /* kernel */
  for(i=0; i<nr; i++) {
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(dfdtheta, 0, j, k)
	  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(dgdtheta, 0, j, k);
	den = 1.0 + xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f, 0, j, k)
	  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g, 0, j, k); /* check the 1/2 factor in front of g */
	scalar3d_set(jacobian, 0, i, j, k, num/den);
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
	  num = 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(dfdtheta, z, j, k)
	    + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(dgdtheta, z, j, k);
	  den =  xi
	    + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, z, j, k)
	    + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g, z, j, k)
	    + beta/alpha;
	  scalar3d_set(jacobian, z, i, j, k, num/den);
	}
      }
    }
  }
  
  /* external domain */
  for(i=0; i<nr; i++) {
    xi = -cos(PI*i/(nr-1));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = 0.25*(xi*xi + xi - 2.0)*scalar2d_get(dfdtheta, nz-1, j, k);
	den = 1.0 + 0.25*(xi*xi + xi - 2.0)*scalar2d_get(f, nz-1, j, k);
	scalar3d_set(jacobian, nz-1, i, j, k, num/den);
      }
    }
  }
  

  /* free memory here */
}


/***************************************************/
/* Take partial derivative with respect to theta.  */
/*           d f / d theta'                        */
/* The prime represents the grid coordinate        */
/* not the physical coordinate system.             */
/***************************************************/
void dfdthetaprime(coeff *dfdt_coeff, coeff *f_coeff)
{
  int imag;
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  
  nz = dfdt_coeff->nz;
  nr = dfdt_coeff->nr;
  nt = dfdt_coeff->nt;
  np = dfdt_coeff->np;
  
  npc = ( np / 2 ) + 1;
  
  /* take d/dtheta derivative of f */
  /* and return them in dfdt */
  
  /* d/dt[a_j cos(j theta)] = -j a_j sin(j theta) */  
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) { 
      for(i=0; i<nr; i++) { 
	for(k=2*imag; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	  for(j=0; j<nt; j++) {
	    coeff_set(dfdt_coeff, z, i, j, k, imag, 
		      -j * coeff_get(f_coeff, z, i, j, k, imag));
	  }
	}
      }
    }
  }
  
  /* d/dt[a_j sin(j theta)] = +j a_j cos(j theta) */
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) {   
      for(i=0; i<nr; i++) { 
	for(k=1; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	  for(j=1; j<nt-1; j++) {  
	    coeff_set(dfdt_coeff, z, i, j, k, imag, 
		      j * coeff_get(f_coeff, z, i, j, k, imag));
	  }
	}
      }
    }
  }
  
}


/*******************************************************************************************/
/* Take f and divide it by sin(theta):                                                     */
/*           1/sin(theta) * f                                                              */
/* This does not allow for the input f to have switched basis functions (thetashift != 0). */
/* This never happens in Jacobian, gradient, or Laplacian so it doesn't matter             */
/*******************************************************************************************/
void dividebysin(coeff *fbysin_coeff, coeff *f_coeff)
{
  int imag;
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  double fbysin_j;
  double f_jp1;
  double fbysin_jp2;
  
  nz = fbysin_coeff->nz;
  nr = fbysin_coeff->nr;
  nt = fbysin_coeff->nt;
  np = fbysin_coeff->np;
  
  npc = ( np / 2 ) + 1;
  
  /* take f_coeff and divide by sin(theta)*/
  /* and return them in fbysin_coeff */
  
  /* divide cos(j theta) series by sin(theta) */  
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) { 
      for(i=0; i<nr; i++) { 
	for(k=2*imag; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	  /* set [1/sin(theta)]f coefficient for cos[(nt-1)theta] to zero: */
	  coeff_set(fbysin_coeff, z, i, nt-1, k, imag, 0.0);
	  /* set [1/sin(theta)]f coefficient for cos[(nt-2)theta]: */
	  fbysin_j = -2.0*coeff_get(f_coeff, z, i, nt-1, k, imag);
	  coeff_set(fbysin_coeff, z, i, nt-2, k, imag, fbysin_j);
	  /* now set [1/sin(theta)]f for j<=nt-3: */
	  for(j=nt-3; j>0; j--) { /* treat j=0 explicitly after loop */
	    f_jp1 = coeff_get(f_coeff, z, i, j+1, k, imag);
	    fbysin_jp2 = coeff_get(fbysin_coeff, z, i, j+2, k, imag);
	    fbysin_j = -2.0*f_jp1 + fbysin_jp2;
	    coeff_set(fbysin_coeff, z, i, j, k, imag, fbysin_j);
	  }
	  coeff_set(fbysin_coeff, z, i, 0, k, imag, 0.0);
	}
      }
    }
  }
  
  /* divide sin(j theta) series by sin(theta) */ 
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) {   
      for(i=0; i<nr; i++) { 
	for(k=1; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	  /* set [1/sin(theta)]f coefficient for sin[(nt-1)theta] to zero: */
	  coeff_set(fbysin_coeff, z, i, nt-1, k, imag, 0.0);
	  /* set [1/sin(theta)]f coefficient for sin[(nt-2)theta]: */
	  fbysin_j = 2.0*coeff_get(f_coeff, z, i, nt-1, k, imag);
	  coeff_set(fbysin_coeff, z, i, nt-2, k, imag, fbysin_j);
	  /* now set [1/sin(theta)]f for j<=nt-3: */
	  for(j=nt-3; j>=0; j--) {
	    f_jp1 = coeff_get(f_coeff, z, i, j+1, k, imag);
	    fbysin_jp2 = coeff_get(fbysin_coeff, z, i, j+2, k, imag);
	    fbysin_j = (2.0*f_jp1 + fbysin_jp2)/(1+delta(j, 0));
	    coeff_set(fbysin_coeff, z, i, j, k, imag, fbysin_j);
	  }
	}
      }
    }
  }
  
}


/*******************************************************************************************/
/* Take f and divide it by sin(theta):                                                     */
/*           1/sin(theta) * f                                                              */
/* This does not allow for the input f to have switched basis functions (thetashift != 0). */
/* This never happens in the Jacobian, gradient, or Laplacian so it doesn't matter.        */
/*******************************************************************************************/
void dividebysin_bound(bound_coeff *fbysin_bound_coeff, bound_coeff *f_bound_coeff)
{
  int imag;
  int z;
  int j;
  int k;
  int nz;
  int nt;
  int np;
  int npc;
  double fbysin_j;
  double f_jp1;
  double fbysin_jp2;
  
  nz = fbysin_bound_coeff->nz;
  nt = fbysin_bound_coeff->nt;
  np = fbysin_bound_coeff->np;
  
  npc = ( np / 2 ) + 1;
  
  /* take f_bound_coeff and divide by sin(theta)*/
  /* and return them in fbysin_bound_coeff */
  
  /* divide cos(j theta) series by sin(theta) */  
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) { 
      for(k=2*imag; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	/* set [1/sin(theta)]f coefficient for cos[(nt-1)theta] to zero: */
	bound_coeff_set(fbysin_bound_coeff, z, nt-1, k, imag, 0.0);
	/* set [1/sin(theta)]f coefficient for cos[(nt-2)theta]: */
	fbysin_j = -2.0*bound_coeff_get(f_bound_coeff, z, nt-1, k, imag);
	bound_coeff_set(fbysin_bound_coeff, z, nt-2, k, imag, fbysin_j);
	/* now set [1/sin(theta)]f for j<=nt-3: */
	for(j=nt-3; j>0; j--) { /* treat j=0 explicitly after loop */
	  f_jp1 = bound_coeff_get(f_bound_coeff, z, j+1, k, imag);
	  fbysin_jp2 = bound_coeff_get(fbysin_bound_coeff, z, j+2, k, imag);
	  fbysin_j = -2.0*f_jp1 + fbysin_jp2;
	  bound_coeff_set(fbysin_bound_coeff, z, j, k, imag, fbysin_j);
	}
	bound_coeff_set(fbysin_bound_coeff, z, 0, k, imag, 0.0);
      }
    }
  }
  
  /* divide sin(j theta) series by sin(theta) */ 
  for(imag=0; imag<=1; imag++) {
    for(z=0; z<nz; z++) {   
      for(k=1; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	/* set [1/sin(theta)]f coefficient for sin[(nt-1)theta] to zero: */
	bound_coeff_set(fbysin_bound_coeff, z, nt-1, k, imag, 0.0);
	/* set [1/sin(theta)]f coefficient for sin[(nt-2)theta]: */
	fbysin_j = 2.0*bound_coeff_get(f_bound_coeff, z, nt-1, k, imag);
	bound_coeff_set(fbysin_bound_coeff, z, nt-2, k, imag, fbysin_j);
	/* now set [1/sin(theta)]f for j<=nt-3: */
	for(j=nt-3; j>=0; j--) {
	  f_jp1 = bound_coeff_get(f_bound_coeff, z, j+1, k, imag);
	  fbysin_jp2 = bound_coeff_get(fbysin_bound_coeff, z, j+2, k, imag);
	  fbysin_j = (2.0*f_jp1 + fbysin_jp2)/(1+delta(j, 0));
	  bound_coeff_set(fbysin_bound_coeff, z, j, k, imag, fbysin_j);
	}
      }
    }
  }
  
}


/**************************************************************/
/* J_3 = 1/Rsin(theta) * dR(xi, theta, phi)/dphi              */
/**************************************************************/
void jacobian3(scalar3d *jacobian, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g)
{
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  double dcoskpdp;
  double dsinkpdp;
  double alpha;
  double beta;
  double xi;
  double num;
  double den;
  bound_coeff *f_fourier;
  bound_coeff *g_fourier;
  bound_coeff *onebysindfdp_fourier;
  bound_coeff *onebysindgdp_fourier;
  scalar2d *onebysindfdp;
  scalar2d *onebysindgdp;
  
  nz = jacobian->nz;
  nr = jacobian->nr;
  nt = jacobian->nt;
  np = jacobian->np;
  
  npc = ( np / 2 ) + 1;
  
  /* allocate memory here */
  f_fourier = bound_coeff_alloc(nz, nt, np);
  g_fourier = bound_coeff_alloc(nz, nt, np);
  onebysindfdp_fourier = bound_coeff_alloc(nz, nt, np);
  onebysindgdp_fourier = bound_coeff_alloc(nz, nt, np);
  onebysindfdp = scalar2d_alloc(nz, nt, np);
  onebysindgdp = scalar2d_alloc(nz, nt, np);
  
  /* go to Fourier space */
  gridtofourier_bound(f_fourier, f);
  gridtofourier_bound(g_fourier, g);
  
  
  /* take d/dphi derivatives of f and g for all zones */
  /* and return them in the same structure. */
  /* Derivatives turn sin to cos and vice versa.  */
  /* d/dp[a_k cos(k phi)] = -k a_k sin(k phi) */
  /* d/dp[a_k sin(k phi)] = +k a_k cos(k phi) */
  for(z=0; z<nz; z++) {
    for(j=0; j<nt; j++) {
      for(k=0; k<npc-1; k++) { /* treat k=npc-1 separately */
	/* f: */
	dcoskpdp = -k * bound_coeff_get(f_fourier, z, j, k, REAL);
	dsinkpdp = k * bound_coeff_get(f_fourier, z, j, k, IMAG);
	bound_coeff_set(f_fourier, z, j, k, REAL, dsinkpdp);
	bound_coeff_set(f_fourier, z, j, k, IMAG, dcoskpdp);
	/* g: */
	dcoskpdp = -k * bound_coeff_get(g_fourier, z, j, k, REAL);
	dsinkpdp = k * bound_coeff_get(g_fourier, z, j, k, IMAG);
	bound_coeff_set(g_fourier, z, j, k, REAL, dsinkpdp);
	bound_coeff_set(g_fourier, z, j, k, IMAG, dcoskpdp);
      }
      bound_coeff_set(f_fourier, z, j, npc-1, REAL, 0.0);
      bound_coeff_set(f_fourier, z, j, npc-1, IMAG, 0.0);
      bound_coeff_set(g_fourier, z, j, npc-1, REAL, 0.0);
      bound_coeff_set(g_fourier, z, j, npc-1, IMAG, 0.0);
    }
  }
  
  /* divide f_fourier and g_fourier by sin(theta) */
  dividebysin_bound(onebysindfdp_fourier, f_fourier);
  dividebysin_bound(onebysindgdp_fourier, g_fourier);

  /* return from Fourier space to grid */
  fouriertogrid_bound(onebysindfdp, onebysindfdp_fourier, 1);
  fouriertogrid_bound(onebysindgdp, onebysindgdp_fourier, 1);
  
  /* kernel */
  for(i=0; i<nr; i++) {
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(onebysindfdp, 0, j, k)
	  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(onebysindgdp, 0, j, k);
	den = 1.0 + xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f, 0, j, k)
	  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g, 0, j, k); /* check the 1/2 factor in front of g */
	scalar3d_set(jacobian, 0, i, j, k, num/den);
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
	  num = 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(onebysindfdp, z, j, k)
	    + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(onebysindgdp, z, j, k);
	  den = xi
	    + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, z, j, k)
	    + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g, z, j, k)
	    + beta/alpha;
	  scalar3d_set(jacobian, z, i, j, k, num/den);
	}
      }
    }
  }
  
  /* external domain */
  for(i=0; i<nr; i++) {
    xi = -cos(PI*i/(nr-1));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = 0.25*(xi*xi+xi-2.0)*scalar2d_get(onebysindfdp, nz-1, j, k);
	den = 1.0 + 0.25*(xi*xi+xi-2.0)*scalar2d_get(f, nz-1, j, k);
	scalar3d_set(jacobian, nz-1, i, j, k, num/den);
      }
    }
  }
  
}


/**************************************************************/
/* 1 / [R sin(theta)] * d^2 R(xi, theta, phi) / [d phi d xi]  */
/**************************************************************/
void onebyrsin_d2rbydpdx(scalar3d *out_grid, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g)
{
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  double dcoskpdp;
  double dsinkpdp;
  double alpha;
  double beta;
  double xi;
  double num;
  double den;
  bound_coeff *f_fourier;
  bound_coeff *g_fourier;
  bound_coeff *onebysindfdp_fourier;
  bound_coeff *onebysindgdp_fourier;
  scalar2d *onebysindfdp;
  scalar2d *onebysindgdp;
  
  nz = out_grid->nz;
  nr = out_grid->nr;
  nt = out_grid->nt;
  np = out_grid->np;
  
  npc = ( np / 2 ) + 1;
  
  /* allocate memory here */
  f_fourier = bound_coeff_alloc(nz, nt, np);
  g_fourier = bound_coeff_alloc(nz, nt, np);
  onebysindfdp_fourier = bound_coeff_alloc(nz, nt, np);
  onebysindgdp_fourier = bound_coeff_alloc(nz, nt, np);
  onebysindfdp = scalar2d_alloc(nz, nt, np);
  onebysindgdp = scalar2d_alloc(nz, nt, np);
  
  /* go to Fourier space */
  gridtofourier_bound(f_fourier, f);
  gridtofourier_bound(g_fourier, g);
  
  
  /* take d/dphi derivatives of f and g for all zones */
  /* and return them in the same structure. */
  /* Derivatives turn sin to cos and vice versa.  */
  /* d/dp[a_k cos(k phi)] = -k a_k sin(k phi) */
  /* d/dp[a_k sin(k phi)] = +k a_k cos(k phi) */
  for(z=0; z<nz; z++) {
    for(j=0; j<nt; j++) {
      for(k=0; k<npc-1; k++) { /* treat k=npc-1 separately */
	/* f: */
	dcoskpdp = -k * bound_coeff_get(f_fourier, z, j, k, REAL);
	dsinkpdp = k * bound_coeff_get(f_fourier, z, j, k, IMAG);
	bound_coeff_set(f_fourier, z, j, k, REAL, dsinkpdp);
	bound_coeff_set(f_fourier, z, j, k, IMAG, dcoskpdp);
	/* g: */
	dcoskpdp = -k * bound_coeff_get(g_fourier, z, j, k, REAL);
	dsinkpdp = k * bound_coeff_get(g_fourier, z, j, k, IMAG);
	bound_coeff_set(g_fourier, z, j, k, REAL, dsinkpdp);
	bound_coeff_set(g_fourier, z, j, k, IMAG, dcoskpdp);
      }
      bound_coeff_set(f_fourier, z, j, npc-1, REAL, 0.0);
      bound_coeff_set(f_fourier, z, j, npc-1, IMAG, 0.0);
      bound_coeff_set(g_fourier, z, j, npc-1, REAL, 0.0);
      bound_coeff_set(g_fourier, z, j, npc-1, IMAG, 0.0);
    }
  }
  
  /* divide f_fourier and g_fourier by sin(theta) */
  dividebysin_bound(onebysindfdp_fourier, f_fourier);
  dividebysin_bound(onebysindgdp_fourier, g_fourier);

  /* return from Fourier space to grid */
  fouriertogrid_bound(onebysindfdp, onebysindfdp_fourier, 1);
  fouriertogrid_bound(onebysindgdp, onebysindgdp_fourier, 1);
  
  /* kernel */
  for(i=0; i<nr; i++) {
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = 12.0*xi*xi*(1.0-xi*xi)*scalar2d_get(onebysindfdp, 0, j, k)
	  + 7.5*xi*(1.0-xi*xi)*scalar2d_get(onebysindgdp, 0, j, k);
	den = 1.0 
	  + xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f, 0, j, k)
	  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g, 0, j, k); /* check the 1/2 factor in front of g */
	scalar3d_set(out_grid, 0, i, j, k, num/den);
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
	  num = 0.75*(xi*xi-1.0)*scalar2d_get(onebysindfdp, z, j, k)
	    + 0.75*(-xi*xi+1.0)*scalar2d_get(onebysindgdp, z, j, k);
	  den = xi
	    + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f, z, j, k)
	    + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g, z, j, k)
	    + beta/alpha;
	  scalar3d_set(out_grid, z, i, j, k, num/den);
	}
      }
    }
  }
  
  /* external domain */
  for(i=0; i<nr; i++) {
    xi = -cos(PI*i/(nr-1));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = 0.75*(xi+1.0)*scalar2d_get(onebysindfdp, nz-1, j, k);
	den = 1.0 + 0.25*(xi*xi+xi-2.0)*scalar2d_get(f, nz-1, j, k);
	scalar3d_set(out_grid, nz-1, i, j, k, num/den);
      }
    }
  }
  
}


/***************************************************/
/* Take partial derivative with respect to phi'    */
/*           d f / d phi'                          */
/* The prime represents the grid coordinate        */
/* not the physical coordinate system.             */
/* The input and output can be the same? VERIFY THIS */
/* But of course the input will then become the output */
/***************************************************/
void dfdphiprime(coeff *dfdp_coeff, coeff *f_coeff)
{
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  double dcoskpdp;
  double dsinkpdp;
  
  nz = dfdp_coeff->nz;
  nr = dfdp_coeff->nr;
  nt = dfdp_coeff->nt;
  np = dfdp_coeff->np;
  
  npc = ( np / 2 ) + 1;
  
  /* take d/dphi derivatives of f and g for all zones */
  /* and return them in the same structure. */
  /* Derivatives turn sin to cos and vice versa.  */
  /* d/dp[a_k cos(k phi)] = -k a_k sin(k phi) */
  /* d/dp[a_k sin(k phi)] = +k a_k cos(k phi) */
  for(z=0; z<nz; z++) {
    for(i=0; i<nr; i++) { 
      for(j=0; j<nt; j++) {
	for(k=0; k<npc-1; k++) {
	  dcoskpdp = -k * coeff_get(f_coeff, z, i, j, k, REAL);
	  dsinkpdp = k * coeff_get(f_coeff, z, i, j, k, IMAG);
	  coeff_set(dfdp_coeff, z, i, j, k, REAL, dsinkpdp);
	  coeff_set(dfdp_coeff, z, i, j, k, IMAG, dcoskpdp);
	}
	/* treat k=npc-1 separately */
	coeff_set(dfdp_coeff, z, i, j, npc-1, REAL, 0.0);
	coeff_set(dfdp_coeff, z, i, j, npc-1, IMAG, 0.0);
      }
    }
  }
}


/*lapf_jp2k = (j+2 > nt-1) ? 0.0 : coeff_get(lapf_coeff, z, i, j+2, k, imag);*/
/**************************************************************************************/
/* Evaluate the angular part of the Laplacian for spherical coordinates:              */
/* Delta_{theta phi} = d^2/dtheta^2 + cot(theta) d/dtheta + 1/sin^2(theta) d^2/dphi^2 */
/**************************************************************************************/
void laplace_ang(coeff *lapf_coeff, coeff *f_coeff)
{
  int imag;
  int z;
  int i;
  int j;
  int k;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  coeff *cott_dfdt_coeff;
  coeff *onebysin2t_d2fdt2_coeff;
  double f_jk; 
  double f_jp2k; /* {j+2, k} coefficient of f */
  double cott_dfdt_jk;
  double cott_dfdt_jp2k;
  double onebysin2t_d2fdt2_jk;
  double onebysin2t_d2fdt2_jp2k;
  double onebysin2t_d2fdt2_jp4k;
  double lapf_jk; /* {j, k} coefficient of \Delta_{\theta\phi}f */
  int jmax;
  int jmin;

  nz = lapf_coeff->nz;
  nr = lapf_coeff->nr;
  nt = lapf_coeff->nt;
  np = lapf_coeff->np;
  
  npc = ( np / 2 ) + 1;
  
  cott_dfdt_coeff = coeff_alloc(nz, nr, nt, np);
  onebysin2t_d2fdt2_coeff = coeff_alloc(nz, nr, nt, np);
  
  /* cos(j theta) part */
  jmax = nt-1;
  jmin = 0;
  for(z=0; z<nz; z++) { 
    for(imag=0; imag<=1; imag++) {
      for(i=0; i<nr; i++) { 
	for(k=2*imag; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	  for(j=jmax; j>=jmin; j--) {
	    f_jk = coeff_get(f_coeff, z, i, j, k, imag);
	    f_jp2k = (j+2 > jmax) ? 0.0 : coeff_get(f_coeff, z, i, j+2, k, imag);
	    cott_dfdt_jp2k = (j+2 > jmax) ? 0.0 : coeff_get(cott_dfdt_coeff, z, i, j+2, k, imag);
	    onebysin2t_d2fdt2_jp2k = (j+2 > jmax) ? 0.0 : coeff_get(onebysin2t_d2fdt2_coeff, z, i, j+2, k, imag);
	    onebysin2t_d2fdt2_jp4k = (j+4 > jmax) ? 0.0 : coeff_get(onebysin2t_d2fdt2_coeff, z, i, j+4, k, imag);
	    /* the recursion relations: */
	    cott_dfdt_jk = (cott_dfdt_jp2k - j*f_jk - (j+2.0)*f_jp2k) 
	      / (1.0+delta(j, 0));
	    onebysin2t_d2fdt2_jk = (4.0*k*k*f_jp2k + 2.0*onebysin2t_d2fdt2_jp2k - onebysin2t_d2fdt2_jp4k) 
	      / (1.0+delta(j, 0));
	    lapf_jk = -j*j*f_jk + cott_dfdt_jk + onebysin2t_d2fdt2_jk;
	    coeff_set(lapf_coeff, z, i, j, k, imag, lapf_jk);
	    coeff_set(cott_dfdt_coeff, z, i, j, k, imag, cott_dfdt_jk);
	    coeff_set(onebysin2t_d2fdt2_coeff, z, i, j, k, imag, onebysin2t_d2fdt2_jk);
	  }
	}
      }
    }
  }
  
  /* sin(j theta) part */ 
  jmax = nt-2;
  jmin = 1;
  for(z=0; z<nz; z++) {   
    for(imag=0; imag<=1; imag++) {
      for(i=0; i<nr; i++) { 
	for(k=1; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	  for(j=jmax; j>=jmin; j--) {
	    f_jk = coeff_get(f_coeff, z, i, j, k, imag);
	    f_jp2k = (j+2 > jmax) ? 0.0 : coeff_get(f_coeff, z, i, j+2, k, imag);
	    cott_dfdt_jp2k = (j+2 > jmax) ? 0.0 : coeff_get(cott_dfdt_coeff, z, i, j+2, k, imag);
	    onebysin2t_d2fdt2_jp2k = (j+2 > jmax) ? 0.0 : coeff_get(onebysin2t_d2fdt2_coeff, z, i, j+2, k, imag);
	    onebysin2t_d2fdt2_jp4k = (j+4 > jmax) ? 0.0 : coeff_get(onebysin2t_d2fdt2_coeff, z, i, j+4, k, imag);
	    /* the recursion relations: */
	    cott_dfdt_jk = (cott_dfdt_jp2k - j*f_jk - (j+2.0)*f_jp2k) 
	      / (1.0+delta(j, 0));
	    onebysin2t_d2fdt2_jk = (4.0*k*k*f_jp2k + 2.0*onebysin2t_d2fdt2_jp2k - onebysin2t_d2fdt2_jp4k) 
	      / (1.0+delta(j, 0));
	    lapf_jk = -j*j*f_jk + cott_dfdt_jk + onebysin2t_d2fdt2_jk;
	    coeff_set(lapf_coeff, z, i, j, k, imag, lapf_jk);
	    coeff_set(cott_dfdt_coeff, z, i, j, k, imag, cott_dfdt_jk);
	    coeff_set(onebysin2t_d2fdt2_coeff, z, i, j, k, imag, onebysin2t_d2fdt2_jk);
	  }
	}
      }
    }
  }
 
  coeff_free(cott_dfdt_coeff);
  coeff_free(onebysin2t_d2fdt2_coeff);
}


/**********************************************************************/
/* Take grid points for each boundary b_zjk and decompose             */
/* them into coefficients of {cos(j theta), sin(j theta)}e^(i k phi). */
/**********************************************************************/
void gridtofourier_bound(bound_coeff *bcoeff, scalar2d *b_scalar2d)
{
  int nz; /* number of zones */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int j, k;
  int imag;
  scalar2d *b_zjk;
  double *in1dphi; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying even theta for fixed phi */
  fftw_plan plan_forward_phi;
  fftw_plan plan_forward_cos;
  fftw_plan plan_forward_sin;
  fftw_complex *out1dphic; /* fft in phi direction for fixed theta */
  double *out1dcos; /* fct in theta direction for fixed phi (k is even) */
  double *out1dsin; /* fst in theta direction for fixed phi (k is odd) */
  
  nz = b_scalar2d->nz;
  nt = b_scalar2d->nt;
  np = b_scalar2d->np;
  
  if(np%2 == 1) 
	printf("np must be even in gridtofourier_bound\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* copy data to new structure so original data is not erased */
  b_zjk = scalar2d_alloc(nz, nt, np);
  scalar2d_memcpy(b_zjk, b_scalar2d);

  /* Set up arrays to hold the data */
  in1dphi = fftw_malloc ( sizeof ( double ) * np );
  in1dcos = fftw_malloc ( sizeof ( double ) * nt );
  in1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  out1dphic = fftw_malloc ( sizeof ( fftw_complex ) * npc );
  out1dcos = fftw_malloc ( sizeof ( double ) * nt );
  out1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );

  /* Create FFTW plans.  The algorithm will be stored in the fftw_plan structure.         */
  /* The arrays can be changed later (must stay same size) but the plan will be the same. */
  /* The plan usually overwrites the data, so set the data after making the plan.         */
  /* It might be faster to figure out how to use FFTW wisdom to speed up planning.        */
  plan_forward_phi = fftw_plan_dft_r2c_1d ( np, in1dphi, out1dphic, FFTW_ESTIMATE );
  plan_forward_cos = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
  plan_forward_sin = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
  
  /* big loop for each boundary */
  for (z=0; z<nz; z++) {
    
    /*>>>>>>>>>>>>>>>>>>>>>>> PHI DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /* transform for phi for fixed theta value (fixed j) */
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < np; k++ ) {
	in1dphi[k] = scalar2d_get(b_zjk, z, j, k);
      }
      fftw_execute ( plan_forward_phi );
      for ( k = 0; k < npc; k++ ) {
	bound_coeff_set( bcoeff, z, j, k, REAL,
			 out1dphic[k][0]*(2.0-delta(k, 0)-delta(k, npc-1))/np );
	bound_coeff_set( bcoeff, z, j, k, IMAG,
			 out1dphic[k][1]*(-1)*(2.0-delta(k, 0)-delta(k, npc-1))/np );
      }
    } 
    
    /*>>>>>>>>>>>>>>>>>>>>>> THETA DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /* cosine transform for even k. */
    for (imag = 0; imag <= 1; imag++) {
      for ( k = 2*imag; k < npc-imag; k += 2 ) { /* k = 0 and k = npc-1 don't have immaginary parts */
	for ( j = 0; j < nt; j++ ) {
	  in1dcos[j] = bound_coeff_get( bcoeff, z, j, k, imag );
	}
	fftw_execute ( plan_forward_cos );
	for ( j = 0; j < nt; j++ ) {
	  bound_coeff_set( bcoeff, z, j, k, imag, 
			   out1dcos[j]*(2.0-delta(j, 0)-delta(j, nt-1))/(2*(nt-1)) );
	}
      } 
    }
    
    /* sine transform for odd k. */
    for (imag = 0; imag <= 1; imag++) {
      for ( k = 1; k < npc-imag; k += 2 ) {
	for ( j = 1; j < nt-1; j++ ) { /* j = 0 and j = nt-1 don't exist for sine transform */
	  in1dsin[j-1] = bound_coeff_get( bcoeff, z, j, k, imag );
	}
	fftw_execute ( plan_forward_sin );
	bound_coeff_set(bcoeff, z, 0, k, imag, 0.0); /* first element is 0 and not included in RODFT00 */
	for ( j = 1; j < nt-1; j++ ) {
	  bound_coeff_set(bcoeff, z, j, k, imag, out1dsin[j-1]/(nt-1));
	}
	bound_coeff_set(bcoeff, z, nt-1, k, imag, 0.0); /* last element is 0 and not included in RODFT00 */
      } 
    }
    
  } /* end of z loop */
  
  /* Delete plan and arrays */
  scalar2d_free(b_zjk);
  fftw_destroy_plan ( plan_forward_phi );
  fftw_destroy_plan ( plan_forward_cos );
  fftw_destroy_plan ( plan_forward_sin );
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
void fouriertogrid_bound(scalar2d *b_zjk, bound_coeff *b_bound_coeff, int thetashift)
{
  int nz; /* number of zones */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int j, k;
  int imag;
  int kstart;
  bound_coeff *bcoeff;
  fftw_complex *in1dphic; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying even theta for fixed phi */
  fftw_plan plan_backward_cos;
  fftw_plan plan_backward_sin;
  fftw_plan plan_backward_phi;
  double *out1dphi; /* inverse fft in phi direction for fixed theta */
  double *out1dcos; /* fct in theta direction for fixed phi (k is even) */
  double *out1dsin; /* fst in theta direction for fixed phi (k is odd) */
  
  nz = b_bound_coeff->nz;
  nt = b_bound_coeff->nt;
  np = b_bound_coeff->np;
  
  if(np%2 == 1)
    printf("np must be even in fouriertogrid_bound\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* copy data to new structure so original data is not erased */
  bcoeff = bound_coeff_alloc(nz, nt, np);
  bound_coeff_memcpy(bcoeff, b_bound_coeff);

  /* Set up arrays to hold the data */
  in1dcos = fftw_malloc ( sizeof ( double ) * nt );
  in1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) ); 
  out1dcos = fftw_malloc ( sizeof ( double ) * nt );
  out1dsin = fftw_malloc ( sizeof ( double ) * (nt-2) );
  in1dphic = fftw_malloc ( sizeof ( fftw_complex ) * npc );
  out1dphi = fftw_malloc ( sizeof ( double ) * np );
  
  /* Create FFTW plans.  The algorithm will be stored in the fftw_plan structure.         */
  /* The arrays can be changed later (must stay same size) but the plan will be the same. */
  /* The plan usually overwrites the data, so set the data after making the plan.         */
  /* It might be faster to figure out how to use FFTW wisdom to speed up planning.        */
  plan_backward_cos = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
  plan_backward_sin = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
  plan_backward_phi = fftw_plan_dft_c2r_1d ( np, in1dphic, out1dphi, FFTW_ESTIMATE );
  
  /* big loop for each boundary */
  for (z=0; z<nz; z++) {
    
    /*>>>>>>>>>>>>>>>>>>>>>> INVERSE THETA TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/

    /* cosine transform */
    /*    
	  kstart = 
	  imag = 0 and thetashift = 0 : 0
	  imag = 1 and thetashift = 0 : 2
	  imag = 0 and thetashift = 1 : 1
	  imag = 1 and thetashift = 1 : 1
    */
    for (imag = 0; imag <= 1; imag++) {
      kstart = (thetashift==0 ? 2*imag : 1);
      for ( k = kstart; k < npc-imag; k += 2 ) { /* start k at 0 for real part and 2 for imag part */
	for ( j = 0; j < nt; j++ ) {
	  in1dcos[j] = bound_coeff_get( bcoeff, z, j, k, imag ) / (2.0-delta(j, 0)-delta(j, nt-1)); /* Correct for my convention */
	}
	fftw_execute ( plan_backward_cos );
	for ( j = 0; j < nt; j++ ) {
	  bound_coeff_set( bcoeff, z, j, k, imag, out1dcos[j] );
	}
      } 
    }
    
    /* sine transform */
    /*    
	  kstart = 
	  imag = 0 and thetashift = 0 : 1
	  imag = 1 and thetashift = 0 : 1
	  imag = 0 and thetashift = 1 : 0
	  imag = 1 and thetashift = 1 : 2
    */
    for (imag = 0; imag <= 1; imag++) {
      kstart = (thetashift==0 ? 1 : 2*imag);
      for ( k = kstart; k < npc-imag; k += 2 ) {
	/* j = 0 and j = nt-1 don't exist for sine transform */
	for ( j = 1; j < nt-1; j++ ) {
	  in1dsin[j-1] = bound_coeff_get( bcoeff, z, j, k, imag ) / 2.0;
	}
	fftw_execute ( plan_backward_sin );
	for ( j = 1; j < nt-1; j++ ) {
	  bound_coeff_set(bcoeff, z, j, k, imag, out1dsin[j-1]);
	}
      } 
    }
    
    /*>>>>>>>>>>>>>>>>>>>>>> INVERSE PHI TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /* transform for phi for fixed theta value (fixed j) */
    for ( j = 0; j < nt; j++ ) {
      for ( k = 0; k < npc; k++ ) {
	in1dphic[k][0] = bound_coeff_get(bcoeff, z, j, k, REAL)
	  / (2.0-delta(k, 0)-delta(k, npc-1));
	in1dphic[k][1] = bound_coeff_get(bcoeff, z, j, k, IMAG)
	  * (-1) / (2.0-delta(k, 0)-delta(k, npc-1));
      }
      fftw_execute ( plan_backward_phi );
      for ( k = 0; k < np; k++ ) {
	scalar2d_set( b_zjk, z, j, k, out1dphi[k] );
      }
    } 
    
  } /* end of z loop */
  
  /* Delete plan and arrays */
  bound_coeff_free(bcoeff);
  fftw_destroy_plan ( plan_backward_cos );
  fftw_destroy_plan ( plan_backward_sin );
  fftw_destroy_plan ( plan_backward_phi );
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
/**********************************************************************/
void gridtofourier(coeff *coeff, scalar3d *c_scalar3d, int xishift, int thetashift)
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int i, j, k;
  int imag;
  scalar3d *c_zijk;
  double *in1dphi; /* picks out varying phi for fixed theta */
  double *in1dcos; /* picks out varying even theta for fixed phi */
  double *in1dsin; /* picks out varying even theta for fixed phi */
  double *in1dxi; /* even part of kernel or all parts of other zones */
  double *in1dxiodd; /* odd part of kernel */
  fftw_plan plan_forward_phi;  
  fftw_plan plan_forward_cos;
  fftw_plan plan_forward_sin;
  fftw_plan plan_forward_xi;
  fftw_plan plan_forward_xiodd;
  fftw_complex *out1dphic; /* fft in phi direction for fixed theta */
  double *out1dcos; /* fct in theta direction for fixed phi (k is even) */
  double *out1dsin; /* fst in theta direction for fixed phi (k is odd) */
  double *out1dxi;
  double *out1dxiodd;
  
  nz = c_scalar3d->nz;
  nr = c_scalar3d->nr;
  nt = c_scalar3d->nt;
  np = c_scalar3d->np;
  
  if(np%2 == 1)
	printf("np must be even in gridtofourier\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* copy data to new structure so original data is not erased */
  c_zijk = scalar3d_alloc(nz, nr, nt, np);
  scalar3d_memcpy(c_zijk, c_scalar3d);

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
  
  /* Create FFTW plans.  The algorithm will be stored in the fftw_plan structure.         */
  /* The arrays can be changed later (must stay same size) but the plan will be the same. */
  /* The plan usually overwrites the data, so set the data after making the plan.         */
  /* It might be faster to figure out how to use FFTW wisdom to speed up planning.        */
  plan_forward_phi = fftw_plan_dft_r2c_1d ( np, in1dphi, out1dphic, FFTW_ESTIMATE );
  plan_forward_cos = fftw_plan_r2r_1d ( nt, in1dcos, out1dcos, FFTW_REDFT00, FFTW_ESTIMATE );
  plan_forward_sin = fftw_plan_r2r_1d ( nt-2, in1dsin, out1dsin, FFTW_RODFT00, FFTW_ESTIMATE );
  plan_forward_xi = fftw_plan_r2r_1d ( nr, in1dxi, out1dxi, FFTW_REDFT00, FFTW_ESTIMATE );
  plan_forward_xiodd = fftw_plan_r2r_1d ( nr-1, in1dxiodd, out1dxiodd, FFTW_REDFT01, FFTW_ESTIMATE ); /* do type-3 DCT (inverse of type-2 DCT) */ 
  
  /* big loop to find coefficients of basis functions for each zone */
  for (z=0; z<nz; z++) {
    
    /*>>>>>>>>>>>>>>>>>>>>>>> PHI DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /* transform for phi for fixed xi and theta */
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  in1dphi[k] = scalar3d_get(c_zijk, z, i, j, k);
	}
	fftw_execute ( plan_forward_phi );
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
    for (imag = 0; imag <= 1; imag++) {
      for ( i = 0; i < nr; i++ ) {
	for ( k = 2*imag; k < npc-imag; k += 2 ) {
	  for ( j = 0; j < nt; j++ ) {
	    in1dcos[j] = coeff_get( coeff, z, i, j, k, imag );
	  }
	  fftw_execute ( plan_forward_cos );
	  for ( j = 0; j < nt; j++ ) {
	    coeff_set( coeff, z, i, j, k, imag, 
		       out1dcos[j]*(2.0-delta(j, 0)-delta(j, nt-1))/(2*(nt-1)) );
	  }
	} 
      }	  
    }
    
    /* sine transform for odd k. */
    for (imag = 0; imag <= 1; imag++) {
      for ( i = 0; i < nr; i++ ) {
	for ( k = 1; k < npc-imag; k += 2 ) {
	  for ( j = 1; j < nt-1; j++ ) { /* j = 0 and j = nt-1 don't exist for sine transform */
	    in1dsin[j-1] = coeff_get( coeff, z, i, j, k, imag );
	  }
	  fftw_execute ( plan_forward_sin );
	  coeff_set(coeff, z, i, 0, k, imag, 0.0); /* first element is 0 and not included in RODFT00 */
	  for ( j = 1; j < nt-1; j++ ) {
	    coeff_set(coeff, z, i, j, k, imag, out1dsin[j-1]/(nt-1));
	  }
	  coeff_set(coeff, z, i, nt-1, k, imag, 0.0); /* last element is 0 and not included in RODFT00 */
	} 
      }
    }
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>> XI DECOMPOSITION <<<<<<<<<<<<<<<<<<<<<<<<<<*/
    /*              transform for xi for fixed theta and phi               */

    if(z == 0) { /* kernel */
      
      /* even Chebyshev series if j is even */
      for (imag = 0; imag <= 1; imag++) {
	for ( k = imag; k < npc-imag; k++ ) { /* k = 0 and k = npc-1 don't have immaginary parts */
	  for(j = xishift; j < nt-xishift; j += 2) {
	    /*for ( j = 0; j < nt; j += 2 ) {*/
	    for ( i = 0; i < nr; i++ ) {
	      in1dxi[i] = coeff_get( coeff, z, i, j, k, imag );
	    }
	    fftw_execute ( plan_forward_xi );
	    for ( i = 0; i < nr; i++ ) {
	      coeff_set( coeff, z, i, j, k, imag, 
			 out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
	    }
	  }
	}
      }
      
      /* odd Chebyshev series if j is odd */ 
      for (imag = 0; imag <= 1; imag++) {
	for ( k = imag; k < npc-imag; k++ ) { /* k = 0 and k = npc-1 don't have immaginary parts */
	  for(j = 1-xishift; j < nt-1+xishift; j += 2) {
	    /*for ( j = 1; j < nt-1; j += 2 ) {*/
	    /* drop 1st point (it's zero) and reverse other points */		
	    for ( i = 1; i < nr; i++ ) {
	      in1dxiodd[i-1] = coeff_get( coeff, z, nr-i, j, k, imag );
	    }
	    /* do type-3 DCT */
	    fftw_execute ( plan_forward_xiodd );
	    for ( i = 0; i < nr-1; i++ ) {
	      coeff_set( coeff, z, i, j, k, imag, 
			 out1dxiodd[i]/(nr-1) );
	    }
	    coeff_set( coeff, z, nr-1, j, k, imag, 0.0 ); /* set coefficient for T_{nr-1}(xi) = 0 */
	  }
	}
      }
      
    } else { /* other zones */
      
      /* normal Chebyshev series */
      for (imag = 0; imag <= 1; imag++) {
	for ( k = imag; k < npc-imag; k++ ) { /* k = 0 and k = npc-1 don't have immaginary parts */
	  for ( j = 0; j < nt; j++ ) {
	    for ( i = 0; i < nr; i++ ) {
	      in1dxi[i] = coeff_get( coeff, z, i, j, k, imag );
	    }
	    fftw_execute ( plan_forward_xi );
	    for ( i = 0; i < nr; i++ ) {
	      coeff_set( coeff, z, i, j, k, imag,
			 out1dxi[i]*neg1toi(i)*(2.0-delta(i, 0)-delta(i, nr-1))/(2*(nr-1)) );
	    }
	  }
	}
      }
      
    } /* end of radial decomposition */
    
  } /* end of z loop */
  
  /* Delete plan and arrays */
  scalar3d_free(c_zijk);
  fftw_destroy_plan ( plan_forward_phi );
  fftw_destroy_plan ( plan_forward_cos );
  fftw_destroy_plan ( plan_forward_sin );
  fftw_destroy_plan ( plan_forward_xi );
  fftw_destroy_plan ( plan_forward_xiodd );
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
/* and evaluate their values on a grid (xi_i, theta_j, phi_k).              */
/****************************************************************************/
void fouriertogrid(scalar3d *c_zijk, coeff *c_coeff, int xishift, int thetashift)
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction. must be even. */
  int npc; /* number of complex numbers in phi direction */
  int z; /* current boundary of zone */
  int imag;
  int i, j, k;
  int kstart;
  coeff *coeff;
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
  
  nz = c_coeff->nz;
  nr = c_coeff->nr;
  nt = c_coeff->nt;
  np = c_coeff->np;
  
  if(np%2 == 1)
	printf("np must be even in gridtofourier\n");
  
  npc = ( np / 2 ) + 1; /* the first and last numbers are real (np re+im values) */
  
  /* copy data to new structure so original data is not erased */
  coeff = coeff_alloc(nz, nr, nt, np);
  coeff_memcpy(coeff, c_coeff);

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
    /*print_coeff(coeff);*/ 
    /*>>>>>>>>>>>>>>>>>>>>>> INVERSE THETA TRANSFORM <<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /* cosine transform */
    /*    
	  kstart = 
	  imag = 0 and thetashift = 0 : 0
	  imag = 1 and thetashift = 0 : 2
	  imag = 0 and thetashift = 1 : 1
	  imag = 1 and thetashift = 1 : 1
    */
    for(imag=0; imag<=1; imag++) {
      kstart = (thetashift==0 ? 2*imag : 1);
      for ( i = 0; i < nr; i++ ) {
	for ( k = kstart; k < npc-imag; k += 2 ) { /* start k at 0 for real part and 2 for imag part */
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
    
    /* sine transform */
    /*    
	  kstart = 
	  imag = 0 and thetashift = 0 : 1
	  imag = 1 and thetashift = 0 : 1
	  imag = 0 and thetashift = 1 : 0
	  imag = 1 and thetashift = 1 : 2
    */
    for(imag=0; imag<=1; imag++) {
      kstart = (thetashift==0 ? 1 : 2*imag);
      for ( i = 0; i < nr; i++ ) {
	for ( k = kstart; k < npc-imag; k += 2 ) {
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
  coeff_free(coeff);
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




/**********************************************************************************/
/****    THE FUNCTIONS BELOW ARE PROBABLY USELESS BUT SAVE THEM ANYWAY.        ****/
/**********************************************************************************/


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
