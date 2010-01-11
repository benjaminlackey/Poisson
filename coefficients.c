/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "poisson.h"


/****************************************************/
/* Allocate memory to store coordinates of 3d grid. */
/****************************************************/
grid3d *grid3d_alloc(int nz, int nr, int nt, int np)
{
  int i;
  grid3d *points;
  double *data;
  
  points = (grid3d *)malloc(sizeof(grid3d));
  
  /* nz*nr*nt*np data points */
  data = (double *)malloc(sizeof(double)*3*nz*nr*nt*np);
  
  for(i=0; i<3*nz*nr*nt*np; i++) {
    data[i] = 0.0;
  }

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
  int i;
  scalar3d *f_zijk;
  double *data;
  
  f_zijk = (scalar3d *)malloc(sizeof(scalar3d));
  
  /* nz*nr*nt*np data points */
  data = (double *)malloc(sizeof(double)*nz*nr*nt*np);

  for(i=0; i<nz*nr*nt*np; i++) {
    data[i] = 0.0;
  }
  
  f_zijk->nz = nz;
  f_zijk->nr = nr;
  f_zijk->nt = nt;
  f_zijk->np = np;
  f_zijk->data = data;
  
  return f_zijk;
}


void scalar3d_free(scalar3d *f)
{
  free(f->data);
  free(f);
}

/**********************************************************/
/* Allocate memory to store data on a 2d grid with zones. */
/**********************************************************/
scalar2d *scalar2d_alloc(int nz, int nt, int np)
{
  int i;
  scalar2d *b_zjk;
  double *data;
  
  b_zjk = (scalar2d *)malloc(sizeof(scalar2d));
  
  /* nz*nt*np data points */
  data = (double *)malloc(sizeof(double)*nz*nt*np);
  
  for(i=0; i<nz*nt*np; i++) {
    data[i] = 0.0;
  }
  
  b_zjk->nz = nz;
  b_zjk->nt = nt;
  b_zjk->np = np;
  b_zjk->data = data;
  
  return b_zjk;
}


/**********************************************************************/
/* Allocate memory to store coefficients in spherical harmonic basis. */
/**********************************************************************/
ylm_coeff *ylm_coeff_alloc(int nz, int nr, int nt, int np)
{
  int i;
  ylm_coeff *c;
  double *data;

  /* sizeof(coeff) is the number of bytes needed to store 4 integers
     and the address of data.  Memory for the actual data is not
     allocated here: */
  c = (ylm_coeff *)malloc(sizeof(ylm_coeff));
  
  /* Memory for the data is allocated here: */
  data = (double *)malloc(sizeof(double)*nz*nr*nt*(nt+1));

  /* initialize the data to 0.0 */
  for(i=0; i<nz*nr*nt*(nt+1); i++) {
    data[i] = 0.0;
  }
  
  c->nz = nz;
  c->nr = nr;
  c->nt = nt;
  c->np = np;
  c->data = data;
  
  return c;
}


void ylm_coeff_free(ylm_coeff *c)
{
  free(c->data);
  free(c);
}


/***********************************************************/
/* Allocate memory to store coefficients in Fourier basis. */
/***********************************************************/
coeff *coeff_alloc(int nz, int nr, int nt, int np)
{
  /* z goes from 0 to nz-1 */
  /* i goes from 0 to nr-1 */
  /* j goes from 0 to nt-1 */
  /* k goes from 0 to np-1 */
 
  int i;
  int npc; /* number of complex numbers in phi direction */
  coeff *c;
  double *data;
  
  npc = np/2 + 1;

  /* sizeof(coeff) is the number of bytes needed to store 4 integers
     and the address of data.  Memory for the actual data is not
     allocated here: */
  c = (coeff *)malloc(sizeof(coeff));

  /* Memory for the data is allocated here: */
  /* There are nz*nr*nt*2*npc doubles. */
  data = (double *)malloc(sizeof(double)*nz*nr*nt*2*npc);
  
  for(i=0; i<nz*nr*nt*2*npc; i++) {
    data[i] = 0.0;
  }
  
  c->nz = nz;
  c->nr = nr;
  c->nt = nt;
  c->np = np;
  c->data = data;
  
  return c;
}


void coeff_free(coeff *c)
{
  free(c->data);
  free(c);
}

/* /\***************************************************\/ */
/* /\* Allocate memory to store boundary coefficients. *\/ */
/* /\***************************************************\/ */
/* bound_coeff *bound_coeff_alloc(int nz, int nl) */
/* { */
/*   /\* z goes from 0 to nz-1 *\/ */
/*   /\* l goes from 0 to lmax=nl-1 *\/ */
/*   /\* m goes from -l to l for l *\/ */
  
/*   bound_coeff *b; */
/*   double *data; */
  
/*   b = (bound_coeff *)malloc(sizeof(bound_coeff)); */
  
/*   data = (double *)malloc(sizeof(double)*nz*nl*nl); */
  
/*   b->nz = nz; */
/*   b->nl = nl; */
/*   b->data = data; */
  
/*   return b; */
/* } */


/**************************************************/
/* Allocate memory to store boundary coefficients */
/* in fourier basis.                              */
/**************************************************/
bound_coeff *bound_coeff_alloc(int nz, int nt, int np)
{
  /* z goes from 0 to nz-1 */
  /* theta goes from 0 to nt-1 */
  /* phi goes from 0 to np-1 */
  
  int i;
  int npc; /* number of complex numbers in phi direction */
  bound_coeff *b;
  double *data;

  npc = np/2 + 1;

  b = (bound_coeff *)malloc(sizeof(bound_coeff));
  
  data = (double *)malloc(sizeof(double)*nz*nt*2*npc);
  
  for(i=0; i<nz*nt*2*npc; i++) {
    data[i] = 0.0;
  }

  b->nz = nz;
  b->nt = nt;
  b->np = np;
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


/**************************************************/
/* Set spherical harmonic coefficient c_ylm to x. */
/**************************************************/
/*extern inline */
void ylm_coeff_set(ylm_coeff *c, int z, int i, int l, int m, int comp, double x)
{
  /* do error checking */
#ifdef DEBUG
  if((z >= c->nz)||(i >= c->nr)||(l >= c->nt)||(m > l)||(m > c->np)||(m < 0)||(comp<0)||(comp>1))
    printf("index out of bounds in ylm_coeff_set\n");
#endif
  
  /* index = z*nr*nt*(nt+1) + i*nt*(nt+1) + l*(l+1) + 2m + comp */
  /*       = (z*nr + i)*nt*(nt+1) + l*(l+1) + 2m + comp */
  c->data[(z*c->nr + i)*c->nt*(c->nt+1) + l*(l+1) + 2*m + comp] = x;
}


/****************************************/
/* Set Fourier coefficient c_zijk to x. */
/****************************************/
/*extern inline */
void coeff_set(coeff *c, int z, int i, int j, int k, int imag, double x)
{ 
  int npc = (c->np)/2 + 1;
  
  /* do error checking */
#ifdef DEBUG
  if((z >= c->nz)||(i >= c->nr)||(j >= c->nt)||(k >= npc)||(imag<0)||(imag>1))
    printf("index out of bounds in coeff_set\n");
#endif
  
  /* index = z*nr*nt*2*npc + i*nt*2*npc + j*2*npc + 2*k + imag */
  /*       = 2(((z*nr + i)*nt + j)*npc + k) + imag */
  c->data[2*(((z*c->nr + i)*c->nt + j)*npc + k) + imag] = x;
}


/* /\********************************************\/ */
/* /\* Set coefficient b_zlm for boundary to x. *\/ */
/* /\********************************************\/ */
/* /\*extern inline *\/ */
/* void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x) */
/* {  */
/*   /\* do error checking *\/ */
/* #ifdef DEBUG */
/*   if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1)) */
/*     printf("index out of bounds in bound_coeff_set\n"); */
/* #endif */
  
/*   /\* index = z*nl^2 + n*nl^2 + L^2 + 2m + i+ (delta_0m - 1) *\/ */
/*   b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1] = x; */
/* } */

/********************************************/
/* Set coefficient b_zjk for boundary to x. */
/********************************************/
/*extern inline */
void bound_coeff_set(bound_coeff *b, int z, int j, int k, int i, double x)
{ 
  int npc = (b->np)/2 + 1;

  /* do error checking */
#ifdef DEBUG
  if((z >= b->nz)||(j >= b->nt)||(k >= npc)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_set\n");
#endif
  
  /* index = 2*z*nt*npc + 2*j*npc + 2*k + i */
  b->data[2*((z*b->nt + j)*npc + k) + i] = x;
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
    printf("index out of bounds in grid3d_get\n");
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
    printf("index out of bounds in scalar3d_get\n");
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
    printf("index out of bounds in scalar2d_get\n");
#endif
  
  /* index = z*nt*np + j*np + k */
  /*       = (z*nt + j)*np + k  */
  return b->data[(z*b->nt + j)*b->np + k];
} 


/********************************************************/
/* Get value of spherical harmonic coefficient c_ylm.   */
/********************************************************/
/*extern inline */
double ylm_coeff_get(ylm_coeff *c, int z, int i, int l, int m, int comp)
{
  /* do error checking */
#ifdef DEBUG
  if((z >= c->nz)||(i >= c->nr)||(l >= c->nt)||(m > l)||(m > c->np)||(m < 0)||(comp<0)||(comp>1))
    printf("index out of bounds in ylm_coeff_get\n");
#endif
  
  /* index = z*nr*nt*(nt+1) + i*nt*(nt+1) + l*(l+1) + 2m + comp */
  /*       = (z*nr + i)*nt*(nt+1) + l*(l+1) + 2m + comp */
  return c->data[(z*c->nr + i)*c->nt*(c->nt+1) + l*(l+1) + 2*m + comp];
}


/*********************************************/
/* Get value of Fourier coefficient c_znlm.  */
/*********************************************/
/*extern inline */
double coeff_get(coeff *c, int z, int i, int j, int k, int imag)
{
  int npc = (c->np)/2 + 1;
  
  /* do error checking */
#ifdef DEBUG
  if((z >= c->nz)||(i >= c->nr)||(j >= c->nt)||(k >= npc)||(imag<0)||(imag>1))
    printf("index out of bounds in coeff_get\n");
#endif  
  
  
  return c->data[2*(((z*c->nr + i)*c->nt + j)*npc + k) + imag];
} 


/* /\*****************************\/ */
/* /\* Get value of b_zlm.      *\/ */
/* /\*****************************\/ */
/* /\*extern inline *\/ */
/* double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i) */
/* { */
/*   /\* do error checking *\/ */
/* #ifdef DEBUG */
/*   if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1)) */
/*     printf("index out of bounds in bound_coeff_get\n"); */
/* #endif */
  
/*   return b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1]; */
/* } */


/*****************************/
/* Get value of b_zlm.      */
/*****************************/
/*extern inline */
double bound_coeff_get(bound_coeff *b, int z, int j, int k, int i)
{ 
  int npc = (b->np)/2 + 1;

  /* do error checking */
#ifdef DEBUG
  if((z >= b->nz)||(j >= b->nt)||(k >= npc)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_get\n");
#endif
  
  return b->data[2*((z*b->nt + j)*npc + k) + i];
}
