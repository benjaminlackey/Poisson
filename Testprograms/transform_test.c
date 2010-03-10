/*********************************************************************************
 *      *
 * Author: Benjamin D. Lackey                                                    *
 *********************************************************************************/

/* To compile type: gcc -I/opt/local/include -I/Users/lackey/Research/Poisson/ -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W /Users/lackey/Research/Poisson/print.c /Users/lackey/Research/Poisson/coefficients.c /Users/lackey/Research/Poisson/coordinatemap.c /Users/lackey/Research/Poisson/fourierylmconversions.c /Users/lackey/Research/Poisson/matrixoperators.c /Users/lackey/Research/Poisson/residual.c /Users/lackey/Research/Poisson/remap.c /Users/lackey/Research/Poisson/gradient.c /Users/lackey/Research/Poisson/radial.c /Users/lackey/Research/Poisson/poisson.h transform_test.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>

/* fftw header */
#include <fftw3.h>

/* own header */
#include "poisson.h"

double boundary(int z, double theta, double phi);
double field(int z, double xi, double theta, double phi);

int main (void)
{
  int nz = 3;
  int nr = 25; /* must be odd? */
  int nt;
  int np = 16; /* must be even */

  scalar2d *boundary_scalar2d;
  scalar2d *f_scalar2d;
  scalar2d *g_scalar2d;
  gsl_vector *alpha_vector;
  gsl_vector *beta_vector;
  scalar3d *field_scalar3d;
 
  bound_coeff *f_bound_coeff;
  bound_coeff *g_bound_coeff;
  coeff *field_coeff;
  
  nt = np/2 + 1;
  
  /* there are nz-1 boundaries but there are still nz of each f and g for each zone */ 
  boundary_scalar2d = scalar2d_alloc(nz-1, nt, np);
  f_scalar2d = scalar2d_alloc(nz, nt, np);
  g_scalar2d = scalar2d_alloc(nz, nt, np);
  alpha_vector = gsl_vector_calloc(nz);
  beta_vector = gsl_vector_calloc(nz);
  field_scalar3d = scalar3d_alloc(nz, nr, nt, np);
  
  f_bound_coeff = bound_coeff_alloc(nz, nt, np);
  g_bound_coeff = bound_coeff_alloc(nz, nt, np);
  field_coeff = coeff_alloc(nz, nr, nt, np);
  
/*   /\* evaluate boundary function on gridpoints *\/ */
/*   boundarytogrid(boundary_scalar2d, boundary); */
  
/*   /\* determine the surface quantities: alpha_vector, beta_vector, f_scalar2d, g_scalar2d *\/ */
/*   map_physicaltogrid(boundary_scalar2d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d); */
/*   print_vector(alpha_vector);  */
/*   print_vector(beta_vector);  */
/* /\*   print_scalar2d(f_scalar2d); *\/ */
/* /\*   print_scalar2d(g_scalar2d); *\/ */
/*   gridtofourier_bound(f_bound_coeff, f_scalar2d); */
/*   gridtofourier_bound(g_bound_coeff, g_scalar2d); */
/*   print_bound_coeff(f_bound_coeff); */
/*   print_bound_coeff(g_bound_coeff); */

/*   /\* evaluate analytical solution at gridpoints *\/ */
/*   functiontogrid(field_scalar3d, alpha_vector, beta_vector, f_scalar2d, g_scalar2d, field); */
/*   gridtofourier(field_coeff, field_scalar3d, 0, 0); */
   /*print_coeff(field_coeff);*/
  
  functiontogrid_xi(field_scalar3d, field);
  gridtofourier(field_coeff, field_scalar3d, 0, 0);
  print_coeff(field_coeff);


  scalar2d_free(boundary_scalar2d);
  scalar2d_free(f_scalar2d);
  scalar2d_free(g_scalar2d);
  gsl_vector_free(alpha_vector);
  gsl_vector_free(beta_vector);
  scalar3d_free(field_scalar3d);

  bound_coeff_free(f_bound_coeff);
  bound_coeff_free(g_bound_coeff);
  coeff_free(field_coeff);

  return 0;
}


double boundary(int z, double theta, double phi)
{
  int L1 = 1;
  int m1 = 0;

  if(z==0)
    return 5.0 + cos(theta);
  else
    return 10.0;
}


double field(int z, double xi, double theta, double phi)
{
  return xi*cos(theta);
}


