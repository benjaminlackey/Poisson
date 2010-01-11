/* To compile type: gcc -g -lm -lgsl -lgslcblas -Wall -pedantic -ansi continuitytest.c */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* constants and simple functions */
#define PI 3.141592653589793
#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
#define cosipiby2(i) ((i)%2 ? 0 : ((i)%4==0 ? 1 : -1)) /* \cos(i\pi/2) */

/* function prototypes */
void solve_continuity(int L, 
					  gsl_vector *boundary, 
					  gsl_matrix *particular_m, 
					  gsl_vector *homo_coeff);

void print_vector(gsl_vector *v);
void print_matrix(gsl_matrix *m);

/* main program */
int main (void)
{
  int N = 5;
  int L = 0;
  
  int nz = 3; /* number of zones (must be >= 2) */

  double b_data[]={5, 10, 20};

  double m_data[]={1, 1, 1, 1, 1,
				   1, 1, 1, 1, 1,
				   0, 0, 0, 0, 0};

  double h_data[]={0, 0, 0, 0};

  gsl_vector_view boundary = gsl_vector_view_array (b_data, nz-1);
  gsl_matrix_view particular_m = gsl_matrix_view_array (m_data, nz, N);
  gsl_vector_view homo_coeff = gsl_vector_view_array (h_data, 2*nz-2);
  
  /*                |Chebyshev coefficients of solution in zone 0   | */
  /*  		        |Chebyshev coefficients of solution in zone 1   | */
  /* particular_m = |                       ...                     | */
  /*                |Chebyshev coefficients of solution in zone nz-1| */
 

  solve_continuity(L, &boundary.vector, &particular_m.matrix, &homo_coeff.vector);


  print_vector(&homo_coeff.vector);
  
  return 0;
}


/********************************************************************/
/* Find coefficients for homogeneous solution to satisfy continuity */
/* given the particular solution.                                   */
/* Returns homo_coeff = (A0, A1, B1, A2, B2, ..., Az-2, Bz-2, Bz-1) */
/********************************************************************/
void solve_continuity(int L, 
					  gsl_vector *boundary, 
					  gsl_matrix *particular_m, 
					  gsl_vector *homo_coeff)
{
  int i;
  int k;
  int z = (*particular_m).size1; /* number of zones (must be >= 2) */
  int N = (*particular_m).size2;
  
  double R;
  double alpha;
  double fin;
  double dfin;
  double fout;
  double dfout;
  
  gsl_matrix *cont_m = gsl_matrix_calloc (2*z-2, 2*z-2); /* set elements to zero */
  gsl_vector *cont_v = gsl_vector_calloc (2*z-2); /* set elements to zero */

  int luint;
  gsl_permutation *permute = gsl_permutation_alloc (2*z-2);
  
  /* z zones */  
  /* z-1 boundaries between zones */
  /* z-2 is index of last boundary in boundary vector */
  /* 0 is index of boundary between kernal and 1st shell */
  /* 1...z-3 are indices of boundaries between two shells */
  /* z-2 is index of boundary between last shell and external domain */
  
  /* set values for cont_m */
  if(z < 2){ /* error */
	
	printf("There must be at least 2 zones");
	
  }else if(z == 2){ /* there is only a kernel and an external domain */
	R = gsl_vector_get(boundary, 0);
	/* row for continuity */
	gsl_matrix_set(cont_m, 0, 0, pow(R, L));
	gsl_matrix_set(cont_m, 0, 1, -pow(R, -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 1, 0, L*pow(R, L-1));
	gsl_matrix_set(cont_m, 1, 1, (L+1)*pow(R, -L-2));
	
  }else{ /* there are shells */
	
	/* internal boundary: */
	R = gsl_vector_get(boundary, 0);
	/* row for continuity */
	gsl_matrix_set(cont_m, 0, 0, pow(R, L));
	gsl_matrix_set(cont_m, 0, 1, -pow(R, L));
	gsl_matrix_set(cont_m, 0, 2, -pow(R, -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 1, 0, L*pow(R, L-1));
	gsl_matrix_set(cont_m, 1, 1, -L*pow(R, L-1));
	gsl_matrix_set(cont_m, 1, 2, (L+1)*pow(R, -L-2));
	
	/* shell boundaries: */
	
	for(i=1; i<=z-3; i++){
	  R = gsl_vector_get(boundary, i);
	  /* row for continuity */
	  gsl_matrix_set(cont_m, 2*i, 2*i-1, pow(R, L));
	  gsl_matrix_set(cont_m, 2*i, 2*i, pow(R, -L-1));
	  gsl_matrix_set(cont_m, 2*i, 2*i+1, -pow(R, L));
	  gsl_matrix_set(cont_m, 2*i, 2*i+2, -pow(R, -L-1));
	  /* row for continuity of 1st derivative */
	  gsl_matrix_set(cont_m, 2*i+1, 2*i-1, L*pow(R, L-1));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i, -(L+1)*pow(R, -L-2));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i+1, -L*pow(R, L-1));
	  gsl_matrix_set(cont_m, 2*i+1, 2*i+2, (L+1)*pow(R, -L-2));	
	}
	
	/* external boundary: */
	R = gsl_vector_get(boundary, z-2);
	/* row for continuity */
	gsl_matrix_set(cont_m, 2*z-4, 2*z-5, pow(R, L));
	gsl_matrix_set(cont_m, 2*z-4, 2*z-4, pow(R, -L-1));
	gsl_matrix_set(cont_m, 2*z-4, 2*z-3, -pow(R, -L-1));
	/* row for continuity of 1st derivative */
	gsl_matrix_set(cont_m, 2*z-3, 2*z-5, L*pow(R, L-1));
	gsl_matrix_set(cont_m, 2*z-3, 2*z-4, -(L+1)*pow(R, -L-2));
	gsl_matrix_set(cont_m, 2*z-3, 2*z-3, (L+1)*pow(R, -L-2));
	
  }
  
  printf("continuity matrix:\n");
  print_matrix(cont_m);
  
  /* set values for cont_v */
  if(z == 2){ /* there is only a kernel and an external domain */
	
	fin = dfin = fout = dfout = 0; 
	R = gsl_vector_get(boundary, 0);
	/* function and derivative on inside: */
	alpha = R;
	if(!(L%2)){ /* L is even */
	  for(k=0; k<N; k++){	
		fin += gsl_matrix_get(particular_m, 0, k);
		dfin += 4*k*k*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	} else { /* L is odd */
	  for(k=0; k<N; k++){
		fin += gsl_matrix_get(particular_m, 0, k);		
		dfin += (2*k+1)*(2*k+1)*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	}
	/* function and derivative on outside: */	
	alpha = -0.5/R;
	for(k=0; k<N; k++){	
	  fout += neg1toi(k)*gsl_matrix_get(particular_m, 1, k);
	  dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, 1, k);
	}
	dfout /= -(alpha*R*R);
	
	gsl_vector_set(cont_v, 0, -fin+fout);
	gsl_vector_set(cont_v, 1, -dfin+dfout);
	
  }else{ /* there are shells */

	/* interface between kernel and first shell */
	fin = dfin = fout = dfout = 0; 
	/* function and derivative on inside: */
	R = gsl_vector_get(boundary, 0);
	alpha = R;
	if(!(L%2)){ /* L is even */
	  for(k=0; k<N; k++){	
		fin += gsl_matrix_get(particular_m, 0, k);
		dfin += 4*k*k*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	} else { /* L is odd */
	  for(k=0; k<N; k++){
		fin += gsl_matrix_get(particular_m, 0, k);		
		dfin += (2*k+1)*(2*k+1)*gsl_matrix_get(particular_m, 0, k);
	  }
	  dfin /= alpha;
	}

	/* function and derivative on outside: */	
	alpha = 0.5*(gsl_vector_get(boundary, 1)-R);
	for(k=0; k<N; k++){	
	  fout += neg1toi(k)*gsl_matrix_get(particular_m, 1, k);
	  dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, 1, k);	
	  printf("k=%d\tfout=%f\tdfout=%f\t(-1)^{k+1}=%d\n", k, fout, dfout, neg1toi(k+1));
	}
	printf("k=%d\tfout=%f\tdfout=%f\t(-1)^{k+1}=%d\n", k, fout, dfout, neg1toi(k+1));
	dfout /= alpha;
	printf("alpha=%f\n", alpha);
	printf("fin=%f\tdfin=%f\tfout=%f\tdfout=%f\n", fin, dfin, fout, dfout);

	gsl_vector_set(cont_v, 0, -fin+fout);
	gsl_vector_set(cont_v, 1, -dfin+dfout);
	
	/* interface between shells */
	/* (skipped if z<=3) */
	for(i=1; i<=z-3; i++){
	  fin = dfin = fout = dfout = 0; 
	  R = gsl_vector_get(boundary, i);
	  /* function and derivative on inside: */
	  alpha = 0.5*(R-gsl_vector_get(boundary, i-1));
	  for(k=0; k<N; k++){	
		fin += gsl_matrix_get(particular_m, i, k);
		dfin += k*k*gsl_matrix_get(particular_m, i, k);
	  }
	  dfin /= alpha;
	  /* function and derivative on outside: */	
	  alpha = 0.5*(gsl_vector_get(boundary, i+1)-R);
	  for(k=0; k<N; k++){	
		fout += neg1toi(k)*gsl_matrix_get(particular_m, i+1, k);
		dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, i+1, k);
	  }
	  dfout /= alpha;
	  
	  gsl_vector_set(cont_v, 2*i, -fin+fout);
	  gsl_vector_set(cont_v, 2*i+1, -dfin+dfout);
	}
	
	/* interface between last shell and external domain */
	fin = dfin = fout = dfout = 0; 
	R = gsl_vector_get(boundary, z-2);
	/* function and derivative on inside: */
	alpha = 0.5*(R-gsl_vector_get(boundary, z-3));
	for(k=0; k<N; k++){	
	  fin += gsl_matrix_get(particular_m, z-2, k);
	  dfin += k*k*gsl_matrix_get(particular_m, z-2, k);
	}
	dfin /= alpha;
	/* function and derivative on outside: */	
	alpha = -0.5/R;
	for(k=0; k<N; k++){	
	  fout += neg1toi(k)*gsl_matrix_get(particular_m, z-1, k);
	  dfout += neg1toi(k+1)*k*k*gsl_matrix_get(particular_m, z-1, k);
	}
	dfout /= -1/(alpha*R*R);
	
	gsl_vector_set(cont_v, 2*z-4, -fin+fout);
	gsl_vector_set(cont_v, 2*z-3, -dfin+dfout);
  }
  
  printf("vector for continuity equation:\n");
  print_vector(cont_v);
  
  /* Solve for homo_coeff in cont_m.homo_coeff = cont_v */
  gsl_linalg_LU_decomp (cont_m, permute, &luint); 
  gsl_linalg_LU_solve (cont_m, permute, cont_v, homo_coeff); 
  
  printf("solution to continuity equation:\n");
  print_vector(homo_coeff);

  /* free memory */
  gsl_matrix_free(cont_m);
  gsl_vector_free(cont_v);
  gsl_permutation_free(permute);
}


/******************************************************/
/*    Functions for printing matrices and vectors.    */
/******************************************************/

void print_vector(gsl_vector *v)
{
  int n;
  int i;

  n = (*v).size;
  
  for (i = 0; i < n; i++){
    printf ("%8g\t", gsl_vector_get (v, i));
  } 
  printf("\n\n");
}


void print_matrix(gsl_matrix *m)
{
  int rows, columns;
  int i, j;

  rows = (*m).size1;
  columns = (*m).size2;
  
  for (i = 0; i < rows; i++){
    for (j = 0; j < columns; j++)
      printf ("%8g\t", gsl_matrix_get (m, i, j));
    printf ("\n");
  } 
  printf("\n");
}
