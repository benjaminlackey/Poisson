/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "poisson.h"

/*******************/
/* Print a vector. */
/*******************/
void print_vector(gsl_vector *v)
{
  int n;
  int i;

  n = (*v).size;
  
  for (i = 0; i < n; i++){
    printf ("%8g  ", gsl_vector_get (v, i));
  } 
  printf("\n\n");
}


/*******************/
/* Print a matrix. */
/*******************/
void print_matrix(gsl_matrix *m)
{
  int rows, columns;
  int i, j;

  rows = (*m).size1;
  columns = (*m).size2;
  
  for (i = 0; i < rows; i++){
    for (j = 0; j < columns; j++) 
	  printf ("%8g  ", gsl_matrix_get (m, i, j));
	/*printf ("%.18e  ", gsl_matrix_get (m, i, j));*/
    printf ("\n");
  } 
  printf("\n");
}

/*****************************/
/* Print all the grid points */
/*****************************/
void print_grid3d(grid3d *points){
  int nz = points->nz;
  int nr = points->nr;
  int nt = points->nt;
  int np = points->np;
  int z, i, j, k;
  
  for(z=0; z<nz; z++) 
	for(i=0; i<nr; i++) 
	  for(j=0; j<nt; j++) 
		for(k=0; k<np; k++)
		  printf("(%d, %f, %f, %f)\t", z, grid3d_get(points, z, i, j, k, 0),
				 grid3d_get(points, z, i, j, k, 1), 
				 grid3d_get(points, z, i, j, k, 2));
  printf("\n");
}


/********************************************************/
/* Print all the values on each surface for a function. */
/********************************************************/
void print_scalar2d(scalar2d *b){
  int nz = b->nz;
  int nt = b->nt;
  int np = b->np;
  int z, j, k;
  
  for(z=0; z<nz; z++) {
    printf("zone %d:\n    ", z);
	for(k=0; k<np; k++)
	  printf("phi_%d:  ", k);
	printf("\n");
	for(j=0; j<nt; j++) {
	  printf("theta_%d:  ", j);
	  for(k=0; k<np; k++)
		printf("%e  ", scalar2d_get(b, z, j, k)); 
	  printf("\n"); 
	}
  }
  printf("\n");
}


/*****************************************************/
/* Print all the values in each zone for a function. */
/*****************************************************/
void print_scalar3d(scalar3d *f){
  int nz = f->nz;
  int nr = f->nr;
  int nt = f->nt;
  int np = f->np;
  int z, i, j, k;
  
  for(z=0; z<nz; z++) {
	for(i=0; i<nr; i++) {
	  printf("z_%d, xi_%d:\n\t", z, i); 
	  for(k=0; k<np; k++)
		printf("phi_%d:\t", k);
	  printf("\n");
	  for(j=0; j<nt; j++) {
		printf("theta_%d:\t", j);
		for(k=0; k<np; k++)
		  printf("%f\t", scalar3d_get(f, z, i, j, k)); 
		printf("\n"); 
	  }
	}
  }
  printf("\n");
}


/*****************************************************************/
/* Print all the spherical harmonic coefficients for a function. */
/*****************************************************************/
void print_ylm_coeff(ylm_coeff *c){
  int nz = c->nz;
  int nr = c->nr;
  int nt = c->nt;
  int z, i, l, m;
  
  for(z=0; z<nz; z++) {
    for(i=0; i<nr; i++) {
      printf("z = %d, i = %d:\n       ", z, i);
      for(m=0; m<nt; m++)
		printf("m = %d                 ", m);
      printf("\n");
      for(l=0; l<nt; l++) {
		printf("l = %d: ", l);
		for(m=0; m<=l; m++){
		  printf("(%f, %f)  ", ylm_coeff_get(c, z, i, l, m, REAL), ylm_coeff_get(c, z, i, l, m, IMAG));
		}
		printf("\n");
      }
	  printf("\n");
    }
  }
  printf("\n");
}

/**********************************************/
/* Print all the coefficients for a function. */
/* Have k, j printed in inner loop.           */
/**********************************************/
void print_coeff(coeff *c){
  int nz = c->nz;
  int nr = c->nr;
  int nt = c->nt;
  int np = c->np;
  int z, i, j, k;
  
  int npc = np/2 + 1;


  for(z=0; z<nz; z++) {
    for(i=0; i<nr; i++) {
      printf("z = %d, i = %d:\n       ", z, i);
      for(k=0; k<npc; k++)
		printf("k = %d                 ", k);
      printf("\n");
      for(j=0; j<nt; j++) {
		printf("j = %d: ", j);
		for(k=0; k<npc; k++){
		  printf("(%.4e, %.4e)  ", coeff_get(c, z, i, j, k, REAL), coeff_get(c, z, i, j, k, IMAG));
		}
		printf("\n");
      }
	  printf("\n");
    }
  }
  printf("\n");
}


/* void print_coeff(coeff *c){ */
/*   int nz = c->nz; */
/*   int nr = c->nr; */
/*   int nt = c->nt; */
/*   int np = c->np; */
/*   int z, i, j, k; */
  
/*   int npc = np/2 + 1; */


/*   for(z=0; z<nz; z++) { */
/*     for(i=0; i<nr; i++) { */
/*       printf("z = %d, i = %d:\n       ", z, i); */
/*       for(k=0; k<npc; k++) */
/* 		printf("k = %d                 ", k); */
/*       printf("\n"); */
/*       for(j=0; j<nt; j++) { */
/* 		printf("j = %d: ", j); */
/* 		for(k=0; k<npc; k++){ */
/* 		  printf("(%f, %f)  ", coeff_get(c, z, i, j, k, REAL), coeff_get(c, z, i, j, k, IMAG)); */
/* 		} */
/* 		printf("\n"); */
/*       } */
/* 	  printf("\n"); */
/*     } */
/*   } */
/*   printf("\n"); */
/* } */

/**********************************************/
/* Print all the coefficients for a function. */
/* Have z, n printed in inner loop.           */
/**********************************************/
void print_coeff_2(coeff *c){
  int nz = c->nz;
  int nr = c->nr;
  int nt = c->nt;
  int np = c->np;
  int z, i, j, k;
  
  int npc = np/2 + 1;

  for(k=0; k<npc; k++) {
	for(j=0; j<nt; j++) {
	  if(!(k%2)) {
		printf("[cos(%d*p), sin(%d*p)]*cos(%d*t):\n     ", k, k, j);
	  } else {
		printf("[cos(%d*p), sin(%d*p)]*sin(%d*t):\n     ", k, k, j);
	  }
      for(z=0; z<nz; z++){
		if((z==0) && (j%2)){
		  for(i=0; i<nr-1; i++)
			printf("T_%d(x)                 ", 2*i+1);
		  printf("\n");
		} else if((z==0) && !(j%2)) {
		  for(i=0; i<nr; i++)
			printf("T_%d(x)                 ", 2*i);
		  printf("\n");
		} else {
		  for(i=0; i<nr; i++)
			printf("T_%d(x)                 ", i);
		  printf("\n");
		}
		printf("z=%d: ", z);
		for(i=0; i<nr; i++) {
		  printf("(%f, %f)  ", coeff_get(c, z, i, j, k, REAL), coeff_get(c, z, i, j, k, IMAG));
		}
		printf("\n");
		if(z==0)
		  printf("     ");
	  }
	  printf("\n");
	}
  }
  printf("\n");
}


/* /\****************************************\/ */
/* /\* Print all the boundary coefficients. *\/ */
/* /\****************************************\/ */
/* void print_bound_coeff(bound_coeff *b){ */
/*   int nz = b->nz; */
/*   int nl = b->nl; */
/*   int z, L, m; */
  
/*   for(z=0; z<nz; z++) { */
/*     printf("zone %d:\n", z); */
/* 	printf("\tm=0\t"); */
/* 	for(m=1; m<nl; m++) */
/* 	  printf("m=%d\t\t", m); */
/* 	printf("\n"); */
/* 	printf("L=0:\t"); */
/* 	printf("(%f, -)\n", bound_coeff_get(b, z, 0, 0, 0)); */
/* 	for(L=1; L<nl; L++) { */
/* 	  printf("L=%d:\t", L); */
/* 	  printf("(%f, -)\t", bound_coeff_get(b, z, L, 0, 0)); */
/* 	  for(m=1; m<=L; m++){ */
/* 		printf("(%f, %f)\t", bound_coeff_get(b, z, L, m, 0), bound_coeff_get(b, z, L, m, 1)); */
/* 	  } */
/* 	  printf("\n"); */
/* 	} */
/*   } */
/*   printf("\n"); */
/* } */


/****************************************/
/* Print all the boundary coefficients. */
/****************************************/
void print_bound_coeff(bound_coeff *b){
  int nz = b->nz;
  int nt = b->nt;
  int np = b->np;
  int npc;
  int z, j, k;
  
  npc = np/2 + 1;
  
  for(z=0; z<nz; z++) {
    printf("zone %d:\n\t", z);
	
	for(k=0; k<npc; k++)
	  printf("k=%d                   ", k);
	
	printf("\n");
	for ( j = 0; j < nt; j++ ) {
	
	  printf("j=%d:\t", j);
	  for ( k = 0; k < npc; k++ ) {
		printf ("(%f, %f)  ", bound_coeff_get(b, z, j, k, 0), bound_coeff_get(b, z, j, k, 1));
	  }
	  printf ( "\n" );	   
	  
	} 
  }
  printf ( "\n" );
}
