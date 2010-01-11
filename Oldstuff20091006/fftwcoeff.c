/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "poisson.h"

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

/* Store coefficients of f(xi, theta, phi) in each zone. */
typedef struct 
{
  int nz;
  int nn;
  int nl;
  double *data;
} coeff;

/* Store boundary coefficients of B(theta, phi) for each boundary. */
typedef struct
{
  int nz;
  int nl;
  double *data;
} bound_coeff;

#define DEBUG

#define PI 3.141592653589793

#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
#define cosipiby2(i) ((i)%2 ? 0 : ((i)%4==0 ? 1 : -1)) /* \cos(i\pi/2) */

/* allocate, set, get */
grid3d *grid3d_alloc(int nz, int nr, int nt, int np);
void grid3d_set(grid3d *points, int z, int i, int j, int k, double xi, double theta, double phi);
double grid3d_get(grid3d *points, int z, int i, int j, int k, int c);

scalar3d *scalar3d_alloc(int nz, int nr, int nt, int np);
void scalar3d_set(scalar3d *f_zijk, int z, int i, int j, int k, double x);
double scalar3d_get(scalar3d *f_zijk, int z, int i, int j, int k);

scalar2d *scalar2d_alloc(int nz, int nt, int np);
void scalar2d_set(scalar2d *b_zjk, int z, int j, int k, double x);
double scalar2d_get(scalar2d *b_zjk, int z, int j, int k);

coeff *coeff_alloc(int z, int n, int L);
void coeff_set(coeff *c, int z, int n, int L, int m, int i, double x);
double coeff_get(coeff *c, int z, int n, int L, int m, int i);

bound_coeff *bound_coeff_alloc(int z, int L);
void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x);
double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i);
/****************************************************/
/* Allocate memory to store coordinates of 3d grid. */
/****************************************************/
grid3d *grid3d_alloc(int nz, int nr, int nt, int np)
{
  grid3d *points;
  double *data;
  
  points = (grid3d *)malloc(sizeof(grid3d));
  
  /* nz*nr*nt*np data points */
  data = (double *)malloc(sizeof(double)*3*nz*nr*nt*np);
  
  points->nz = nz;
  points->nr = nr;
  points->nt = nt;
  points->np = np;
  points->data = data;
  
  return points;
}

/**********************************************************/
/* Allocate memory to store data on a 3d grid with zones. */
/**********************************************************/
scalar3d *scalar3d_alloc(int nz, int nr, int nt, int np)
{
  scalar3d *f_zijk;
  double *data;
  
  f_zijk = (scalar3d *)malloc(sizeof(scalar3d));
  
  /* nz*nr*nt*np data points */
  data = (double *)malloc(sizeof(double)*nz*nr*nt*np);
  
  f_zijk->nz = nz;
  f_zijk->nr = nr;
  f_zijk->nt = nt;
  f_zijk->np = np;
  f_zijk->data = data;
  
  return f_zijk;
}

/**********************************************************/
/* Allocate memory to store data on a 2d grid with zones. */
/**********************************************************/
scalar2d *scalar2d_alloc(int nz, int nt, int np)
{
  scalar2d *b_zjk;
  double *data;
  
  b_zjk = (scalar2d *)malloc(sizeof(scalar2d));
  
  /* nz*nr*nt*np data points */
  data = (double *)malloc(sizeof(double)*nz*nt*np);
  
  b_zjk->nz = nz;
  b_zjk->nt = nt;
  b_zjk->np = np;
  b_zjk->data = data;
  
  return b_zjk;
}

/*****************************************/
/* Allocate memory to store coefficients */
/* and return a pointer to the memory.   */
/*****************************************/
coeff *coeff_alloc(int nz, int nn, int nl)
{
  /* z goes from 0 to nz-1 */
  /* n goes from 0 to nn-1 */
  /* l goes from 0 to lmax=nl-1 */
  /* m goes from 0 to l */
  /* there are 2*l+1 values for fixed z, n, l */
  
  coeff *c;
  double *data;
  
  /* sizeof(coeff) is the number of bytes needed to store 3 integers
     and the address of data.  Memory for the actual data is not
     allocated here: */
  c = (coeff *)malloc(sizeof(coeff));
  
  /* Memory for the data is allocated here: */
  /* If -l <= m <= l for each l, there are (lmax+1)^2 elements
	 for each n and z (or nl^2 elements where nl=lmax+1). */
  data = (double *)malloc(sizeof(double)*nz*nn*nl*nl);
  
  c->nz = nz;
  c->nn = nn;
  c->nl = nl;
  c->data = data;
  
  return c;
}


/**************************************************/
/* Allocate memory to store boundary coefficients */
/* and return a pointer to the memory.            */
/**************************************************/
bound_coeff *bound_coeff_alloc(int nz, int nl)
{
  /* z goes from 0 to nz-1 */
  /* l goes from 0 to lmax=nl-1 */
  /* m goes from -l to l for l */
  
  bound_coeff *b;
  double *data;
  
  b = (bound_coeff *)malloc(sizeof(bound_coeff));
  
  data = (double *)malloc(sizeof(double)*nz*nl*nl);
  
  b->nz = nz;
  b->nl = nl;
  b->data = data;
  
  return b;
}

/**********************************/
/* Set a coordinate on a 3d grid. */
/**********************************/
/*extern inline */
void grid3d_set(grid3d *points, int z, int i, int j, int k, double xi, double theta, double phi)
{ 
  /* do error checking */
#ifdef DEBUG
  if((z >= points->nz)||(i >= points->nr)||(j >= points->nt)||(k >= points->np))
    printf("index out of bounds in grid3d_set\n");
#endif
  
  /* index = 3*(z*nr*nt*np + i*nt*np + j*np + k) + (0, 1, 2) */
  /*       = 3*(((z*nr + i)*nt + j)*np + k) + (0, 1, 2)      */
  points->data[3*(((z*points->nr + i)*points->nt + j)*points->np + k)] = xi;
  points->data[3*(((z*points->nr + i)*points->nt + j)*points->np + k)+1] = theta; 
  points->data[3*(((z*points->nr + i)*points->nt + j)*points->np + k)+2] = phi;
}

/*******************************/
/* Set value at 3d grid point. */
/*******************************/
/*extern inline */
void scalar3d_set(scalar3d *f, int z, int i, int j, int k, double x)
{ 
  /* do error checking */  
#ifdef DEBUG
  if((z >= f->nz)||(i >= f->nr)||(j >= f->nt)||(k >= f->np))
    printf("index out of bounds in scalar3d_set\n");
#endif
  
  /* index = z*nr*nt*np + i*nt*np + j*np + k */
  /*       = ((z*nr + i)*nt + j)*np + k      */
  f->data[((z*f->nr + i)*f->nt + j)*f->np + k] = x;
}

/*******************************/
/* Set value at 2d grid point. */
/*******************************/
/*extern inline */
void scalar2d_set(scalar2d *b, int z, int j, int k, double x)
{ 
  /* do error checking */
#ifdef DEBUG
  if((z >= b->nz)||(j >= b->nt)||(k >= b->np))
    printf("index out of bounds in scalar2d_set\n");
#endif
  
  /* index = z*nt*np + j*np + k */
  /*       = (z*nt + j)*np + k  */
  b->data[(z*b->nt + j)*b->np + k] = x;
}


/********************************/
/* Set coefficient c_znlm to x. */
/********************************/
/*extern inline */
void coeff_set(coeff *c, int z, int n, int L, int m, int i, double x)
{ 
  /* do error checking */
#ifdef DEBUG
  if((z >= c->nz)||(n >= c->nn)||(L >= c->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in coeff_set\n");
#endif
  
  /* index = z*nn*nl^2 + n*nl^2 + L^2 + 2m + i + (delta_{0m} - 1) */
  /*       = ((z*nn + n)*nl^2 + L^2 + 2m + i + (delta_{0m} - 1) */
  c->data[(z*c->nn + n)*c->nl*c->nl + L*L + 2*m + i + delta(0, m) - 1] = x;
}


/********************************************/
/* Set coefficient b_zlm for boundary to x. */
/********************************************/
/*extern inline */
void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x)
{ 
  /* do error checking */
#ifdef DEBUG
  if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_set\n");
#endif
  
  /* index = z*nl^2 + n*nl^2 + L^2 + 2m + i+ (delta_0m - 1) */
  b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1] = x;
}

/****************************************************/
/* Get value of component of coordinate on 3d grid. */
/****************************************************/
/*extern inline */
double grid3d_get(grid3d *points, int z, int i, int j, int k, int c)
{
  /* do error checking */
#ifdef DEBUG
  if((z >= points->nz)||(i >= points->nr)||(j >= points->nt)
	 ||(k >= points->np)||(c >= 3))
    printf("index out of bounds in coeff_set\n");
#endif
  
  return points->data[3*(((z*points->nr + i)*points->nt + j)*points->np + k) + c];
}
 
/**********************************/
/* Get value of point on 3d grid. */
/**********************************/
/*extern inline */
double scalar3d_get(scalar3d *f, int z, int i, int j, int k)
{
  /* do error checking */
#ifdef DEBUG
  if((z >= f->nz)||(i >= f->nr)||(j >= f->nt)||(k >= f->np))
    printf("index out of bounds in coeff_set\n");
#endif
  
  return f->data[((z*f->nr + i)*f->nt + j)*f->np + k];
} 

/**********************************/
/* Get value of point on 2d grid. */
/**********************************/
/*extern inline */
double scalar2d_get(scalar2d *b, int z, int j, int k)
{
  /* do error checking */
#ifdef DEBUG
  if((z >= b->nz)||(j >= b->nt)||(k >= b->np))
    printf("index out of bounds in coeff_set\n");
#endif
  
  /* index = z*nt*np + j*np + k */
  /*       = (z*nt + j)*np + k  */
  return b->data[(z*b->nt + j)*b->np + k];
} 

/*****************************/
/* Get value of c_znlm.      */
/*****************************/
/*extern inline */
double coeff_get(coeff *c, int z, int n, int L, int m, int i)
{
  /* do error checking */
#ifdef DEBUG
  if((z >= c->nz)||(n >= c->nn)||(L >= c->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in coeff_get\n");
#endif
  
  return c->data[(z*c->nn + n)*c->nl*c->nl + L*L + 2*m + i + delta(0, m) - 1];
} 

/*****************************/
/* Get value of b_zlm.      */
/*****************************/
/*extern inline */
double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i)
{
  /* do error checking */
#ifdef DEBUG
  if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_get\n");
#endif
  
  return b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1];
}
