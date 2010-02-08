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

/******************************************************************************/
/* Evaluate the grotesque quantity that made 6 months of your life dissapear. */
/******************************************************************************/
void residual(scalar3d *residual_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k;
  int nz, nr, nt, np;
  
  scalar3d *j1_scalar3d;
  scalar3d *j2_scalar3d;
  scalar3d *j3_scalar3d;
  scalar3d *rby_xiplusbbya_scalar3d;
  scalar3d *df_dxi_scalar3d;
  scalar3d *onebyr_df_dxi_scalar3d;
  scalar3d *d2r_dxi2_scalar3d;
  scalar3d *onebyr_d2f_dthetadxi_scalar3d;
  scalar3d *onebyr_d2r_dthetadxi_scalar3d;
  scalar3d *onebyrsintheta_d2f_dphidxi_scalar3d;
  scalar3d *onebyrsintheta_d2r_dphidxi_scalar3d;
  scalar3d *onebyrsq_anglaplacef_scalar3d;
  scalar3d *onebyrsq_anglaplacer_scalar3d;
  
  double j1, j2, j3;
  double rby_xiplusbbya;
  double df_dxi;
  double onebyr_df_dxi_d;
  double d2r_dxi2;
  double onebyr_d2f_dthetadxi_d;
  double onebyr_d2r_dthetadxi_d;
  double onebyrsintheta_d2f_dphidxi_d;
  double onebyrsintheta_d2r_dphidxi_d;
  double onebyrsq_anglaplacef_d;
  double onebyrsq_anglaplacer_d;
  
  double rby_xiplusbbya_sq;
  double one_plus_j2sq_plus_j3sq;
  double term1, term2, term3a, term3b, term3;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;
  
  j1_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  j3_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  rby_xiplusbbya_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  df_dxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_df_dxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  d2r_dxi2_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_d2f_dthetadxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyr_d2r_dthetadxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsintheta_d2f_dphidxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsintheta_d2r_dphidxi_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsq_anglaplacef_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  onebyrsq_anglaplacer_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* evaluate each term */
  jacobian1(j1_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  jacobian2(j2_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  jacobian3(j3_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  r_xiplusb_a(rby_xiplusbbya_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  dfdxi_ongrid(df_dxi_scalar3d, f_scalar3d);
  onebyr_df_dxi(onebyr_df_dxi_scalar3d, f_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  d2rdxi2(d2r_dxi2_scalar3d, alpha_vector, f_scalar2d, g_scalar2d);
  onebyr_d2f_dthetadxi(onebyr_d2f_dthetadxi_scalar3d, f_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyr_d2r_dthetadxi(onebyr_d2r_dthetadxi_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); /* you had used f instead of r */
  onebyrsintheta_d2f_dphidxi(onebyrsintheta_d2f_dphidxi_scalar3d, f_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyrsintheta_d2r_dphidxi(onebyrsintheta_d2r_dphidxi_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyrsq_anglaplacef(onebyrsq_anglaplacef_scalar3d, f_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  onebyrsq_anglaplacer(onebyrsq_anglaplacer_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* combine the terms */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  j1 = scalar3d_get(j1_scalar3d, z, i, j, k); /* J_1 */
	  j2 = scalar3d_get(j2_scalar3d, z, i, j, k); /* J_2 */
	  j3 = scalar3d_get(j3_scalar3d, z, i, j, k); /* J_3 */
	  
	  rby_xiplusbbya = scalar3d_get(rby_xiplusbbya_scalar3d, z, i, j, k); /* R/(xi + beta/alpha) */
	  df_dxi = scalar3d_get(df_dxi_scalar3d, z, i, j, k); /* df/dxi */
	  onebyr_df_dxi_d = scalar3d_get(onebyr_df_dxi_scalar3d, z, i, j, k); /* (1/R) * df/dxi */
	  d2r_dxi2 = scalar3d_get(d2r_dxi2_scalar3d, z, i, j, k); /* d^2R/dxi^2 */	 
	  
	  onebyr_d2f_dthetadxi_d = scalar3d_get(onebyr_d2f_dthetadxi_scalar3d, z, i, j, k); /* (1/R) * d^2f/(dtheta dxi) */
	  onebyr_d2r_dthetadxi_d = scalar3d_get(onebyr_d2r_dthetadxi_scalar3d, z, i, j, k); /* (1/R) * d^2R/(dtheta dxi) */
	  
	  onebyrsintheta_d2f_dphidxi_d = scalar3d_get(onebyrsintheta_d2f_dphidxi_scalar3d, z, i, j, k); /* 1/(R sin(theta)) * d^2f/(dphi dxi) */
	  onebyrsintheta_d2r_dphidxi_d = scalar3d_get(onebyrsintheta_d2r_dphidxi_scalar3d, z, i, j, k); /* 1/(R sin(theta)) * d^2R/(dphi dxi) */
	  
	  onebyrsq_anglaplacef_d = scalar3d_get(onebyrsq_anglaplacef_scalar3d, z, i, j, k); /* (1/R^2) * nabla_{theta phi}f */
	  onebyrsq_anglaplacer_d = scalar3d_get(onebyrsq_anglaplacer_scalar3d, z, i, j, k); /* (1/R^2) * nabla_{theta phi}R */
	  
	  rby_xiplusbbya_sq = rby_xiplusbbya*rby_xiplusbbya;
	  one_plus_j2sq_plus_j3sq = 1.0 + j2*j2 + j3*j3;
	  
	  term1 = (rby_xiplusbbya*one_plus_j2sq_plus_j3sq/j1 - 1.0) * 2.0*onebyr_df_dxi_d/j1;
	  
	  term2 = (rby_xiplusbbya_sq*one_plus_j2sq_plus_j3sq/(j1*j1) - 1.0) * onebyrsq_anglaplacef_d;
	  
	  term3a = 2.0*(j2*onebyr_d2f_dthetadxi_d + j3*onebyrsintheta_d2f_dphidxi_d);
	  
	  term3b = (onebyrsq_anglaplacer_d 
		    + (d2r_dxi2*one_plus_j2sq_plus_j3sq/j1 
		       - 2.0*(j2*onebyr_d2r_dthetadxi_d + j3*onebyrsintheta_d2r_dphidxi_d)) / j1) * df_dxi;
	  
	  term3 = (term3a + term3b)/j1;
	  
	  scalar3d_set(residual_scalar3d, z, i, j, k, term1 + term2 + term3);
	  /* printf("z=%d, i=%d, j=%d, k=%d, residual = %.18e\n", z, i, j, k, term1 + term2 + term3); */
	  /* printf("z=%d, i=%d, j=%d, k=%d, term1=%.18e, term2=%.18e, term3a=%.18e, term3b=%.18e\n", z, i, j, k, term1, term2, term3a, term3b); */
	  /* printf("z=%d, i=%d, j=%d, k=%d, r/(xi+b/a)=%.18e\n", z, i, j, k, rby_xiplusbbya); */
	  /* printf("z=%d, i=%d, j=%d, k=%d, (1/R^2)anglapR=%.18e, d^2R/dxi^2=%.18e, (1/R)d^2R/dtdx=%.18e, (1/Rsint)d^2R/dpdx=%.18e\n", z, i, j, k,
	              onebyrsq_anglaplacer_d, d2r_dxi2, onebyr_d2r_dthetadxi_d, onebyrsintheta_d2r_dphidxi_d); */
	}
      }
    }
  }
  
  scalar3d_free(j1_scalar3d);
  scalar3d_free(j2_scalar3d);
  scalar3d_free(j3_scalar3d);
  scalar3d_free(rby_xiplusbbya_scalar3d);
  scalar3d_free(df_dxi_scalar3d);
  scalar3d_free(onebyr_df_dxi_scalar3d);
  scalar3d_free(d2r_dxi2_scalar3d);
  scalar3d_free(onebyr_d2f_dthetadxi_scalar3d);
  scalar3d_free(onebyr_d2r_dthetadxi_scalar3d);
  scalar3d_free(onebyrsintheta_d2f_dphidxi_scalar3d);
  scalar3d_free(onebyrsintheta_d2r_dphidxi_scalar3d);
  scalar3d_free(onebyrsq_anglaplacef_scalar3d);
  scalar3d_free(onebyrsq_anglaplacer_scalar3d);
}


void dfdxi_ongrid(scalar3d *df_dxi_scalar3d, scalar3d *f_scalar3d)
{
  int nz, nr, nt, np;
  coeff *f_coeff;
  coeff *df_dxi_coeff;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;
  
  f_coeff = coeff_alloc(nz, nr, nt, np);
  df_dxi_coeff = coeff_alloc(nz, nr, nt, np);
  
  /* evaluate coefficients of f */
  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  
  /* df/dxi. xishift: 0->1 */
  dfdxi(df_dxi_coeff, f_coeff);
  
  /* go back to gridpoints */
  fouriertogrid(df_dxi_scalar3d, df_dxi_coeff, 1, 0);
  
  coeff_free(f_coeff);
  coeff_free(df_dxi_coeff);
}


void onebyr_df_dxi(scalar3d *onebyr_df_dxi_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int nz, nr, nt, np;
  coeff *f_coeff;
  coeff *df_dxi_coeff;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;
  
  f_coeff = coeff_alloc(nz, nr, nt, np);
  df_dxi_coeff = coeff_alloc(nz, nr, nt, np);
  
  /* evaluate coefficients of f */
  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  
  /* df/dxi. xishift: 0->1 */
  dfdxi(df_dxi_coeff, f_coeff);
  
  /* (1/R) df/dxi */
  dividebyr(onebyr_df_dxi_scalar3d, df_dxi_coeff, 1, 0, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  coeff_free(f_coeff);
  coeff_free(df_dxi_coeff);
}


void onebyrsq_anglaplacef(scalar3d *lapfbyr2_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  int z, i, j, k, imag;
  int nz, nr, nt, np, npc;
  coeff *f_coeff;
  coeff *lapf_coeff;
  bound_coeff *lapfbyxi2_bound_coeff;
  scalar2d *lapfbyxi2_scalar2d;
  scalar3d *r_scalar3d;
  double sum;
  double alpha;
  double roru;
  
  nz = f_scalar3d->nz;
  nr = f_scalar3d->nr;
  nt = f_scalar3d->nt;
  np = f_scalar3d->np;
  
  npc = ( np / 2 ) + 1;
  
  f_coeff = coeff_alloc(nz, nr, nt, np);
  lapf_coeff = coeff_alloc(nz, nr, nt, np);
  lapfbyxi2_bound_coeff = bound_coeff_alloc(1, nt, np);
  lapfbyxi2_scalar2d = scalar2d_alloc(1, nt, np);
  r_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  /* evaluate coefficients of f */
  gridtofourier(f_coeff, f_scalar3d, 0, 0);
  
  /* evaluate angular part of Laplacian in spherical coordinates */
  laplace_ang(lapf_coeff, f_coeff);
  
  /*>>>>>>>>>>> EVALUATE anglapF/xi^2 AT xi=R=0 <<<<<<<<<<<<<<*/
  
  /* Divide by xi^2 at xi = 0. */
  /* Use l'Hopital's rule twice to take lim_{xi->0} lapf(xi, theta, phi)/xi^2. */
  /* Only even Chebyshev polynomials contribute. */
  z = 0;
  for(imag=0; imag<=1; imag++) {
    for(k = imag; k < npc-imag; k++) { /* imaginary parts for k=0 and k=npc-1 are zero */ 
      for(j=0; j < nt; j += 2) {
	sum = 0.0;
	for(i=nr-1; i>=0; i--) {
	  sum += neg1toi(i+1)*4*i*i*coeff_get(lapf_coeff, z, i, j, k, imag);
	}
	bound_coeff_set(lapfbyxi2_bound_coeff, z, j, k, imag, sum/2.0);
      }
    }
  }
  
  /* evaluate at the gridpoints */
  fouriertogrid_bound(lapfbyxi2_scalar2d, lapfbyxi2_bound_coeff, 0);
  /* print_scalar2d(fbyxi_scalar2d); */ 
  
  /*>>>>>>>>>>>> EVALUATE anglapF/R^2 EVERYWHERE <<<<<<<<<<<<<*/
  
  /* convert to values on grid: */
  /* lapfbyr2_grid is still just lapf at grid points */
  fouriertogrid(lapfbyr2_scalar3d, f_coeff, 0, 0);
  
  /* determine radius of each gridpoint */
  rofxtp(r_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d);
  
  /* now evaluate at every gridpoint */
  for ( z = 0; z < nz; z++ ) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nt; j++ ) {
	for ( k = 0; k < np; k++ ) {
	  if(z==0 && i==0) { /* R = 0 */
	    alpha = gsl_vector_get(alpha_vector, z);
	    scalar3d_set(lapfbyr2_scalar3d, z, i, j, k, 
			 scalar2d_get(lapfbyxi2_scalar2d, z, j, k) / (alpha*alpha));
	  } else if (z==nz-1 && i==nr-1) { /* R = inf (U = 0) */
	    scalar3d_set(lapfbyr2_scalar3d, z, i, j, k, 0.0);
	  } else {
	    roru = scalar3d_get(r_scalar3d, z, i, j, k);
	    scalar3d_set(lapfbyr2_scalar3d, z, i, j, k, 
			 scalar3d_get(lapfbyr2_scalar3d, z, i, j, k) / (roru*roru));
	  }
	}
      }
    }
  }
  
  coeff_free(f_coeff);
  coeff_free(lapf_coeff);
  bound_coeff_free(lapfbyxi2_bound_coeff);
  scalar2d_free(lapfbyxi2_scalar2d);
  scalar3d_free(r_scalar3d);
}


void onebyrsq_anglaplacer(scalar3d *out_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d)
{
  double z, i, j, k;
  double nz, nr, nt, np;
  bound_coeff *f_bound_coeff;
  bound_coeff *g_bound_coeff;
  bound_coeff *lapf_bound_coeff;
  bound_coeff *lapg_bound_coeff;
  scalar2d *lapf_scalar2d;
  scalar2d *lapg_scalar2d;
  double xi;
  double alpha, beta;
  double num, den;
  
  nz = out_scalar3d->nz;
  nr = out_scalar3d->nr;
  nt = out_scalar3d->nt;
  np = out_scalar3d->np;
  
  f_bound_coeff = bound_coeff_alloc(nz, nt, np);
  g_bound_coeff = bound_coeff_alloc(nz, nt, np);
  lapf_bound_coeff = bound_coeff_alloc(nz, nt, np);
  lapg_bound_coeff = bound_coeff_alloc(nz, nt, np);
  lapf_scalar2d = scalar2d_alloc(nz, nt, np);
  lapg_scalar2d = scalar2d_alloc(nz, nt, np);
  
  /* evaluate coefficients of f and g */
  gridtofourier_bound(f_bound_coeff, f_scalar2d);
  gridtofourier_bound(g_bound_coeff, g_scalar2d);
  
  /* evaluate angular part of Laplacian in spherical coordinates */
  laplace_ang_bound(lapf_bound_coeff, f_bound_coeff);
  laplace_ang_bound(lapg_bound_coeff, g_bound_coeff);
  
  /* go back to gridpoints */
  fouriertogrid_bound(lapf_scalar2d, lapf_bound_coeff, 0);
  fouriertogrid_bound(lapg_scalar2d, lapg_bound_coeff, 0);
  
  /* kernel */
  z = 0;
  alpha = gsl_vector_get(alpha_vector, z);
  for(i=0; i<nr; i++) {
    xi = sin(PI*i/(2*(nr-1)));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(lapf_scalar2d, z, j, k)
	  + 0.5*xi*(5.0-3.0*xi*xi)*scalar2d_get(lapg_scalar2d, z, j, k);
	den = 1.0 
	  + xi*xi*xi*(3.0-2.0*xi*xi)*scalar2d_get(f_scalar2d, z, j, k)
	  + 0.5*xi*xi*(5.0-3.0*xi*xi)*scalar2d_get(g_scalar2d, z, j, k);
	scalar3d_set(out_scalar3d, z, i, j, k, num / (alpha*den*den));
      }
    }
  }
  
  /* shells */
  for(z=1; z<nz-1; z++) {
    alpha = gsl_vector_get(alpha_vector, z);
    beta = gsl_vector_get(beta_vector, z);
    for(i=0; i<nr; i++) {
      xi = -cos(PI*i/(nr-1));
      for(j=0; j<nt; j++) {
	for(k=0; k<np; k++) {
	  num = alpha*(0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(lapf_scalar2d, z, j, k)
		       + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(lapg_scalar2d, z, j, k));
	  den = alpha*(xi
		       + 0.25*(xi*xi*xi-3.0*xi+2.0)*scalar2d_get(f_scalar2d, z, j, k)
		       + 0.25*(-xi*xi*xi+3.0*xi+2.0)*scalar2d_get(g_scalar2d, z, j, k))
	    + beta;
	  scalar3d_set(out_scalar3d, z, i, j, k, num / (den*den));
	}
      }
    }
  }
  
  /* external domain */
  z = nz-1;
  alpha = gsl_vector_get(alpha_vector, z);
  for(i=0; i<nr; i++) {
    xi = -cos(PI*i/(nr-1));
    for(j=0; j<nt; j++) {
      for(k=0; k<np; k++) {
	num = 0.25*(xi+2.0)*scalar2d_get(lapf_scalar2d, z, j, k);
	den = 1.0 + 0.25*(xi*xi+xi-2.0)*scalar2d_get(f_scalar2d, z, j, k);
	scalar3d_set(out_scalar3d, z, i, j, k, num / (alpha*den*den));
      }
    }
  }
  
  bound_coeff_free(f_bound_coeff);
  bound_coeff_free(g_bound_coeff);
  bound_coeff_free(lapf_bound_coeff);
  bound_coeff_free(lapg_bound_coeff);
  scalar2d_free(lapf_scalar2d);
  scalar2d_free(lapg_scalar2d);
}


/*****************************************************************************************/
/* Evaluate the angular part of the Laplacian for spherical coordinates on the boundary: */
/* Delta_{theta phi} = d^2/dtheta^2 + cot(theta) d/dtheta + 1/sin^2(theta) d^2/dphi^2    */
/*****************************************************************************************/
void laplace_ang_bound(bound_coeff *lapf_bound_coeff, bound_coeff *f_bound_coeff)
{
  int z, j, k, imag;
  int nz, nt, np, npc;
  bound_coeff *cott_dfdt_bound_coeff;
  bound_coeff *onebysin2t_d2fdt2_bound_coeff;
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
  
  nz = lapf_bound_coeff->nz;
  nt = lapf_bound_coeff->nt;
  np = lapf_bound_coeff->np;
  
  npc = ( np / 2 ) + 1;
  
  cott_dfdt_bound_coeff = bound_coeff_alloc(nz, nt, np);
  onebysin2t_d2fdt2_bound_coeff = bound_coeff_alloc(nz, nt, np);
  
  /* cos(j theta) part */
  jmax = nt-1;
  jmin = 0;
  for(z=0; z<nz; z++) { 
    for(imag=0; imag<=1; imag++) {
      for(k=2*imag; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j=jmax; j>=jmin; j--) {
	  f_jk = bound_coeff_get(f_bound_coeff, z, j, k, imag);
	  f_jp2k = (j+2 > jmax) ? 0.0 : bound_coeff_get(f_bound_coeff, z, j+2, k, imag);
	  cott_dfdt_jp2k = (j+2 > jmax) ? 0.0 : bound_coeff_get(cott_dfdt_bound_coeff, z, j+2, k, imag);
	  onebysin2t_d2fdt2_jp2k = (j+2 > jmax) ? 0.0 : bound_coeff_get(onebysin2t_d2fdt2_bound_coeff, z, j+2, k, imag);
	  onebysin2t_d2fdt2_jp4k = (j+4 > jmax) ? 0.0 : bound_coeff_get(onebysin2t_d2fdt2_bound_coeff, z, j+4, k, imag);
	  /* the recursion relations: */
	  cott_dfdt_jk = (cott_dfdt_jp2k - j*f_jk - (j+2.0)*f_jp2k) 
	    / (1.0+delta(j, 0));
	  onebysin2t_d2fdt2_jk = (4.0*k*k*f_jp2k + 2.0*onebysin2t_d2fdt2_jp2k - onebysin2t_d2fdt2_jp4k) 
	    / (1.0+delta(j, 0));
	  lapf_jk = -j*j*f_jk + cott_dfdt_jk + onebysin2t_d2fdt2_jk;
	  bound_coeff_set(lapf_bound_coeff, z, j, k, imag, lapf_jk);
	  bound_coeff_set(cott_dfdt_bound_coeff, z, j, k, imag, cott_dfdt_jk);
	  bound_coeff_set(onebysin2t_d2fdt2_bound_coeff, z, j, k, imag, onebysin2t_d2fdt2_jk);
	}
      }
    }
  }
  
  /* sin(j theta) part */ 
  jmax = nt-2;
  jmin = 1;
  for(z=0; z<nz; z++) {   
    for(imag=0; imag<=1; imag++) {
      for(k=1; k<npc-imag; k+=2) { /* imaginary parts for k=0 and k=npc-1 are zero */
	for(j=jmax; j>=jmin; j--) {
	  f_jk = bound_coeff_get(f_bound_coeff, z, j, k, imag);
	  f_jp2k = (j+2 > jmax) ? 0.0 : bound_coeff_get(f_bound_coeff, z, j+2, k, imag);
	  cott_dfdt_jp2k = (j+2 > jmax) ? 0.0 : bound_coeff_get(cott_dfdt_bound_coeff, z, j+2, k, imag);
	  onebysin2t_d2fdt2_jp2k = (j+2 > jmax) ? 0.0 : bound_coeff_get(onebysin2t_d2fdt2_bound_coeff, z, j+2, k, imag);
	  onebysin2t_d2fdt2_jp4k = (j+4 > jmax) ? 0.0 : bound_coeff_get(onebysin2t_d2fdt2_bound_coeff, z, j+4, k, imag);
	  /* the recursion relations: */
	  cott_dfdt_jk = (cott_dfdt_jp2k - j*f_jk - (j+2.0)*f_jp2k) 
	    / (1.0+delta(j, 0));
	  onebysin2t_d2fdt2_jk = (4.0*k*k*f_jp2k + 2.0*onebysin2t_d2fdt2_jp2k - onebysin2t_d2fdt2_jp4k) 
	    / (1.0+delta(j, 0));
	  lapf_jk = -j*j*f_jk + cott_dfdt_jk + onebysin2t_d2fdt2_jk;
	  bound_coeff_set(lapf_bound_coeff, z, j, k, imag, lapf_jk);
	  bound_coeff_set(cott_dfdt_bound_coeff, z, j, k, imag, cott_dfdt_jk);
	  bound_coeff_set(onebysin2t_d2fdt2_bound_coeff, z, j, k, imag, onebysin2t_d2fdt2_jk);
	}
      }
    }
  }
  
  bound_coeff_free(cott_dfdt_bound_coeff);
  bound_coeff_free(onebysin2t_d2fdt2_bound_coeff);
}
