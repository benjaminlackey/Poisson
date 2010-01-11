#include <gsl/gsl_matrix.h>
#include <fftw3.h>

/* Store coordinates of grid points */
typedef struct 
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction per zone */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction */
  double *data;
} grid3d;

/* Store data at grid points */
typedef struct 
{
  int nz; /* number of zones */
  int nr; /* number of points in radial direction per zone */
  int nt; /* number of points in theta direction */
  int np; /* number of points in phi direction */
  double *data;
} scalar3d;

/*  */
typedef struct 
{
  int nz;
  int nt;
  int np;
  double *data;
} scalar2d;

/* Store coefficients of f(xi, theta, phi) in each zone */
/* in spherical harmonic basis. */
typedef struct 
{
  int nz;
  int nr;
  int nt;
  int np;
  double *data;
} ylm_coeff;

/* Store coefficients of f(xi, theta, phi) in each zone */
/* in spherical harmonic basis. */
typedef struct 
{
  int nz;
  int nr;
  int nt;
  int np;
  double *data;
} coeff;

/* Store boundary coefficients of B_z(theta, phi) for each boundary */
/* in Fourier basis. */
typedef struct
{
  int nz;
  int nt;
  int np;
  double *data;
} bound_coeff;

#define DEBUG

#define REAL 0 /* Real part of a scalar, coefficient, etc. */
#define IMAG 1
#define PI 3.141592653589793
#define G 6.67259e-8
#define C 2.99792458e10

#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
#define cosipiby2(i) ((i)%2 ? 0 : ((i)%4==0 ? 1 : -1)) /* \cos(i\pi/2) */
#define ABS(i) ((i)<0 ? -(i) : (i))


grid3d *grid3d_alloc(int nz, int nr, int nt, int np);
void grid3d_set(grid3d *points, int z, int i, int j, int k, double xi, double theta, double phi);
double grid3d_get(grid3d *points, int z, int i, int j, int k, int c);


scalar3d *scalar3d_alloc(int nz, int nr, int nt, int np);
void scalar3d_free(scalar3d *f);
void scalar3d_set(scalar3d *f_zijk, int z, int i, int j, int k, double x);
double scalar3d_get(scalar3d *f_zijk, int z, int i, int j, int k);
void scalar3d_addconstant(scalar3d *in_grid, double number, scalar3d *sum_grid);
void scalar3d_add(scalar3d *in1_grid, scalar3d *in2_grid, scalar3d *sum_grid);
void scalar3d_multiply(scalar3d *in1_grid, scalar3d *in2_grid, scalar3d *product_grid);


scalar2d *scalar2d_alloc(int nz, int nt, int np);
void scalar2d_set(scalar2d *b_zjk, int z, int j, int k, double x);
double scalar2d_get(scalar2d *b_zjk, int z, int j, int k);


ylm_coeff *ylm_coeff_alloc(int nz, int nr, int nt, int np);
void ylm_coeff_free(ylm_coeff *c);
void ylm_coeff_set(ylm_coeff *c, int z, int i, int l, int m, int comp, double x);
double ylm_coeff_get(ylm_coeff *c, int z, int i, int l, int m, int comp);


coeff *coeff_alloc(int nz, int nr, int nt, int np);
void coeff_free(coeff *c);
void coeff_set(coeff *c, int z, int i, int j, int k, int imag, double x);
double coeff_get(coeff *c, int z, int i, int j, int k, int imag);


bound_coeff *bound_coeff_alloc(int nz, int nt, int np);
void bound_coeff_set(bound_coeff *b, int z, int j, int k, int i, double x);
double bound_coeff_get(bound_coeff *b, int z, int j, int k, int i);


/* matrix operators */
void set_x(gsl_matrix *m);
void set_xinv(gsl_matrix *m);
void set_dbydx(gsl_matrix *m);
void set_xmin1inv(gsl_matrix *m);

void set_A_kernel_even(gsl_matrix *A_even, int L);
void set_A_kernel_odd(gsl_matrix *m, int L);
void set_A_shell(gsl_matrix *m, int L, double alpha, double beta);
void set_A_ext(gsl_matrix *m, int L);


/* print vectors, matrices, coefficients, etc. */
void print_vector(gsl_vector *v);
void print_matrix(gsl_matrix *m);
void print_grid3d(grid3d *points);
void print_scalar2d(scalar2d *b);
void print_scalar3d(scalar3d *f);
void print_ylm_coeff(ylm_coeff *c);
void print_coeff(coeff *c);
void print_coeff_2(coeff *c);
void print_bound_coeff(bound_coeff *b);


/* functions for mapping spherical coordinates to surface matched coordinates */
void map_physicaltogrid(scalar2d *b_zjk, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f_grid, scalar2d *g_grid);
void map_physicaltogrid_kernel(scalar2d *b_zjk, gsl_vector *alphalist, scalar2d *fodd_grid, scalar2d *geven_grid);
void map_physicaltogrid_shell(scalar2d *b_zjk, int z, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *fin_grid, scalar2d *gout_grid);
void map_physicaltogrid_ext(scalar2d *b_zjk, gsl_vector *alphalist, scalar2d *f_grid);
void rofxtp(scalar3d *rgrid, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g);
void functiontogrid(scalar3d *func_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d, double (*func)(int z, double r, double theta, double phi));

/* fast transforms between gridpoints and basis function coefficients (for grids and boundaries) */
void gridtofourier_bound(bound_coeff *bcoeff_zjk, scalar2d *b_zjk);
void fouriertogrid_bound(scalar2d *b_zjk, bound_coeff *bcoeff, int thetashift);
void gridtofourier(coeff *coeff, scalar3d *c_zijk, int xishift, int thetashift);
void fouriertogrid(scalar3d *c_zijk, coeff *coeff, int xishift, int thetashift);

/* jacobians etc. for converting derivatives between coordinate systems */
void dividebyr(scalar3d *fbyr_grid, coeff *f_coeff, int xishift, int thetashift, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g);
void jacobian1(scalar3d *jacobian, gsl_vector *alphalist, scalar2d *f, scalar2d *g);
void d2rdxi2(scalar3d *rxx, gsl_vector *alphalist, scalar2d *f, scalar2d *g);
void dfdxi(coeff *dfdxi, coeff *f);
void jacobian2(scalar3d *jacobian, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g);
void dfdthetaprime(coeff *dfdt_coeff, coeff *f_coeff);
void dividebysin(coeff *fbysin_coeff, coeff *f_coeff);
void dividebysin_bound(bound_coeff *fbysin_bound_coeff, bound_coeff *f_bound_coeff);
void jacobian3(scalar3d *jacobian, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g);
void onebyrsin_d2rbydpdx(scalar3d *out_grid, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f, scalar2d *g);
void dfdphiprime(coeff *dfdp_coeff, coeff *f_coeff);
void laplace_ang(coeff *lapf_coeff, coeff *f_coeff);

/* calculate the gradient of a scalar */
void gradient_r(scalar3d *gradf_r_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d);
void gradient_theta(scalar3d *gradf_theta_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d);
void gradient_phi(scalar3d *gradf_phi_scalar3d, scalar3d *f_scalar3d, gsl_vector *alpha_vector, gsl_vector *beta_vector, scalar2d *f_scalar2d, scalar2d *g_scalar2d);

/* evaluate integrals used in fourier <--> spherical harmonic conversions */
double int_pmm_cosjt_sint(int m, int j);
double int_pmm_sinjt_sint(int m, int j);
double int_pmm_cosjt(int m, int j);
double int_pmm_sinjt(int m, int j);

/* make individual matrices for fourier <--> spherical harmonic conversions for fixed m */
void set_fouriertoylm_meven(gsl_matrix *transform, int m, int nt);
void set_fouriertoylm_modd(gsl_matrix *transform, int m, int nt);
void set_ylmtofourier_meven(gsl_matrix *transform, int m, int nt);
void set_ylmtofourier_modd(gsl_matrix *transform, int m, int nt);

/* allocate memory for matrices, set their values, free memory */
/* for fourier <--> spherical harmonic conversions */
gsl_matrix **fouriertoylm_matrix_alloc(int nt);
void fouriertoylm_matrix_set(gsl_matrix **fouriertoylm);
void fouriertoylm_matrix_free(gsl_matrix **fouriertoylm);
gsl_matrix **ylmtofourier_matrix_alloc(int nt);
void ylmtofourier_matrix_set(gsl_matrix **ylmtofourier);
void ylmtofourier_matrix_free(gsl_matrix **ylmtofourier);

/* perform fourier <--> spherical harmonic conversions for all m */
void transform_fouriertoylm(coeff *c_fourier, ylm_coeff *c_ylm, gsl_matrix **fouriertoylm);
void transform_ylmtofourier(ylm_coeff *c_ylm, coeff *c_fourier, gsl_matrix **ylmtofourier);


/* functions for making grid, decomposing into basis functions */
void legendre_zero_weight(gsl_vector *zeros, gsl_vector *weights);
/*double source(int zone, double r, double theta, double phi);
  double boundary(int zone, double theta, double phi);*/
double p_lm(int l, int m, double x);
double ytheta_lm(int l, int m, double x);

void makegrid3d(grid3d *x_zijk);
void boundtogrid(scalar2d *b_zjk, grid3d *points, double (*bound)(int, double, double));
void ftogrid(scalar3d *f_zijk, scalar2d *b_zjk, grid3d *points, double (*func)(int, double, double, double));
/*void decompose_bound(bound_coeff *b_zlm, scalar2d *b_zjk);*/
void decompose(coeff *c_znlm, scalar3d *f_zijk);

/* allocate memory for matrices in radial part of solution, set their values, free memory */
gsl_matrix ***radial_matrix_alloc(int nz, int nr, int nt);
void radial_matrix_set(int nz, int nt, gsl_vector *alpha_v, gsl_vector *beta_v, gsl_matrix ***radial);
void radial_matrix_free(int nz, int nt, gsl_matrix ***radial);

/* solve radial part of Poisson equation */
void solve_kernel_even(int L, gsl_matrix *A_even_source, double alpha, gsl_vector *s_even, gsl_vector *f_even);
void solve_kernel_odd(int L, gsl_matrix *A_odd_source, double alpha, gsl_vector *s_odd, gsl_vector *f_odd);
void solve_shell(int L, gsl_matrix *A_shell_source, double alpha, double beta, gsl_vector *s_old, gsl_vector *f);
void solve_ext(int L, gsl_matrix *A_ext_source, double alpha, gsl_vector *s_old, gsl_vector *f);
void solve_radial_particular(int L, gsl_matrix ***radial_matrix, gsl_vector *alpha_vector, gsl_vector *beta_vector, gsl_matrix *source_matrix, gsl_matrix *particular_matrix);
void solve_radial_homogeneous(int L, gsl_vector *alpha_v, gsl_vector *beta_v, gsl_matrix *particular_m, gsl_vector *homo_grow_v, gsl_vector *homo_decay_v);
void homogeneoustochebyshev(ylm_coeff *homo_grow_ylm_coeff, ylm_coeff *homo_decay_ylm_coeff, gsl_vector *alpha_vector, gsl_vector *beta_vector, ylm_coeff *homogeneous_ylm_coeff);
void solve_poisson_spherical(ylm_coeff *field_ylm_coeff, ylm_coeff *source_ylm_coeff, gsl_vector *alpha_v, gsl_vector *beta_v);

/* remap surface and function to new coordinate system */
void findnewsurface(scalar2d *newsurface_scalar2d, scalar3d *enthalpy_scalar3d, gsl_vector *alphalist, gsl_vector *betalist, scalar2d *f_scalar2d, scalar2d *g_scalar2d);
double chebyshevinterpolation(gsl_vector *coeff_vector, double xi);
double spectralinterpolation(coeff *func_coeff, int z, double xi, double theta, double phi);
double roruofxi(int nz, int z, double alpha, double beta, double f, double g, double xi);
double xiofroru(int nz, int z, double alpha, double beta, double f, double g, double r);
void remapgrid(scalar3d *funcnew_scalar3d, gsl_vector *alphanew_vector, gsl_vector *betanew_vector, scalar2d *fnew_scalar2d, scalar2d *gnew_scalar2d, scalar3d *funcold_scalar3d, gsl_vector *alphaold_vector, gsl_vector *betaold_vector, scalar2d *fold_scalar2d, scalar2d *gold_scalar2d);
