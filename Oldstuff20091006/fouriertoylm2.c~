/* To compile type: gcc -lm -lfftw3 -lgsl -lgslcblas fouriertoylm.c print.c coefficients.c poisson.h */

/* c headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "poisson.h"

#define PI 3.141592653589793

#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
#define ABS(i) ((i)<0 ? -(i) : (i))

double int_pmm_cosjt_sint(int m, int j);
double int_cost_pmp1m_cosjt_sint(int m, int j);

int main (void)
{
  int l, m, j;
  double plmint, plminus1mint, plminus2mint;
  int nt = 9;
  
  gsl_matrix *Imatrix;
  gsl_matrix *Jmatrix;
  
  
  for( m = 0; m <= 10; m += 2) {
	for(j = 0; j <= 10; j += 2) {
	  printf("%e  ", int_pmm_cosjt_sint(m, j));
	}
	printf("\n");
  }
  
  printf("\n");
  
  for( m = 0; m <= 10; m += 2) {
	for(j = 0; j <= 10; j += 2) {
	  printf("%e  ", int_cost_pmp1m_cosjt_sint(m, j));
	}
	printf("\n");
  }
  
  printf("\n");
  
  printf("%d, %d\n", 2*ABS(-2-1), ABS(2));
  
  /*  m = 4; */
/*   j = 0; */
/*   plminus2mint = int_pmm_cosjt_sint(m, j); */
/*   printf("%e  ", plminus2mint); */
/*   costplminus1mint = int_cost_pmp1m_cosjt_sint(m, j); */
/*   printf("%e  ", plminus1mint); */
/*   for( l = m+2; l < nt; l++ ) { */
/* 	plmint = (2.0*l-1.0)/(l-m)*plminus1mint - (l+m-1.0)/(l-m)*plminus2mint; */
/* 	printf("%e  ", plmint); */
/* 	plminus2mint = plminus1mint; */
/* 	plminus1mint = plmint; */
/*   } */
/*   printf("\n"); */

  m = 4;
  Imatrix = gsl_matrix_alloc(nt, 2*nt-m);
  Jmatrix = gsl_matrix_alloc(nt+1, 2*nt-m);
  gsl_matrix_set_zero(Imatrix);
  gsl_matrix_set_zero(Jmatrix);

  /* fill matrix for even values of j */
  
  /* fill mth row of I */
  for( j = 0; j < 2*nt-m; j += 2) {
	gsl_matrix_set(Imatrix, m, j, int_pmm_cosjt_sint(m, j));
  }

  /* fill m+1 row of complimentary matrix J */
  for( j = 0; j < 2*nt-m; j += 2) {
	gsl_matrix_set(Jmatrix, m+1, j, int_cost_pmp1m_cosjt_sint(m, j));
  }

  printf("hello I am here\n");
  print_matrix(Imatrix);
  print_matrix(Jmatrix);

  /* fill the other even rows of I and odd rows of j */
  for( l = m+2; l < nt; l += 2) {
	for( j = 0; j < 2*nt-m-2; j += 2) {
	  
	  /* reflect indices of matrix if needed */
	  /* j = 0-2 goes to j = 2 */

	  gsl_matrix_set(Imatrix, l, j, 
					 ((2.0*l-1.0)*gsl_matrix_get(Jmatrix, l-1, j) 
					  - (l+m-1.0)*gsl_matrix_get(Imatrix, l-2, j))
					 /(l-m)
					 );	 
	}
	
	for( j = 0; j < 2*nt-m-2; j += 2) {
	  gsl_matrix_set(Jmatrix, l+1, j, 
					 (2.0*l+1.0)/(l+1.0-m)
					 *(0.5*gsl_matrix_get(Imatrix, l, j)
					   + 0.25*gsl_matrix_get(Imatrix, l, ABS(j-2))
					   + 0.25*gsl_matrix_get(Imatrix, l, j+2))
					 
					 - (l+m)/(l+1.0-m)*gsl_matrix_get(Jmatrix, l-1, j)
					 );
	  
	  printf("j=%d, l+1=%d, %f, %f, %f, %f  \n", j, l+1, 
			 (2.0*l+1.0)/(l+1.0-m), 
			 gsl_matrix_get(Imatrix, l, j), 
			 gsl_matrix_get(Imatrix, l, ABS(j-2)), 
			 gsl_matrix_get(Imatrix, l, j+2)
			 );
	}
	printf("\n");
  }
  
  print_matrix(Imatrix);
  print_matrix(Jmatrix);
  
/*   /\* fill the odd rows of I and even rows of j *\/ */
/*   for( l = m+2; l < nt; l += 2) { */
/* 	for( j = 0; j < 2*nt-m-2; j += 2) { */
/* 	  /\* reflect indices of matrix if needed *\/ */
/* 	  /\* j = 0-2 goes to j = 2 *\/ */

/* 	  gsl_matrix_set(Imatrix, l, j,  */
/* 					 ((2.0*l-1.0)*gsl_matrix_get(Jmatrix, l-1, j)  */
/* 					  - (l+m-1.0)*gsl_matrix_get(Imatrix, l-2, j)) */
/* 					 /(l-m) */
/* 					 );	  */
/* 	} */
	
/* 	for( j = 0; j < 2*nt-m-2; j += 2) { */
/* 	  gsl_matrix_set(Jmatrix, l+1, j,  */
/* 					 (2.0*l+1.0)/(l+1.0-m) */
/* 					 *(0.5*gsl_matrix_get(Imatrix, l, j) */
/* 					   + 0.25*gsl_matrix_get(Imatrix, l, ABS(j-2)) */
/* 					   + 0.25*gsl_matrix_get(Imatrix, l, j+2)) */
					 
/* 					 - (l+m)/(l+1.0-m)*gsl_matrix_get(Jmatrix, l-1, j) */
/* 					 ); */
	  
/* 	  printf("j=%d, l+1=%d, %f, %f, %f, %f  \n", j, l+1,  */
/* 			 (2.0*l+1.0)/(l+1.0-m),  */
/* 			 gsl_matrix_get(Imatrix, l, j),  */
/* 			 gsl_matrix_get(Imatrix, l, ABS(j-2)),  */
/* 			 gsl_matrix_get(Imatrix, l, j+2) */
/* 			 ); */
/* 	} */
/* 	printf("\n"); */
/*   } */
 
  
  return 0;
}


/****************************************************/
/*                                                  */
/*                                                  */
/*                                                  */
/****************************************************/
/* void set_fouriertoylm_meven_jeven(gsl_matrix *matrix) */
/* { */
/*   int n; /\*nXn matrix*\/ */
/*   int i;  */

/*   n = matrix->size1; */

/*   /\*initialize matrix*\/ */
/*   gsl_matrix_set_zero(m); */

/*   /\* pmmint = \int_0^\pi P_m^m(\cos(\theta))\cos(j\theta) d\theta *\/ */
/*   if(j < m) { */
	
/* 	pmmint = neg1toi(j/2)*2; */
/* 	for(i = m+j+3; i <= 2m-1; i += 2) { /\* (2m-1)!!/(m+j+1)!! *\/ */
/*       pmmint *= i; */
/* 	}	 */
/* 	for(i = m-j+2; i <= m+1; i++) { /\* (m+1)!/(m-j+1)!! *\/ */
/*       pmmint *= i; */
/* 	}	 */

/*   } else if(j = m) { */

/* 	pmmint = neg1toi(j/2)*2; */
/* 	for(i = 2; i <= m; i += 2) { /\* m!! *\/ */
/*       pmmint *= i; */
/* 	} */
/* 	fact = 1; */
/* 	for(i = m+3; i <= 2*m+1; i += 2) { /\* (m+3)(m+5)...(2m+1) *\/ */
/*       fact *= i; */
/* 	} */
/* 	pmmint /= fact; */

/*   } else if(j = m+2) { */
	
/* 	pmmint = neg1toi(j/2)*2 / ((2*m+1)*(2*m+3)); */
/*     for(i = 2; i <= m+1; i++) { /\* (m+1)! *\/ */
/*       pmmint *= i; */
/* 	}	 */

/*   }	else { */
	
/* 	pmmint = neg1toi(j-1-m/2)*2; */
/* 	for(i = 2; i <= m+1; i++) { /\* (m+1)! *\/ */
/*       pmmint *= i; */
/* 	} */
/* 	fact = 1; */
/* 	for(i = 2*m+1; i <= m+j+1; i += 2) { /\* (2m+1)(2m+3)...(m+j+1) *\/ */
/*       fact *= i; */
/* 	} */
/* 	pmmint /= fact; */
/* 	fact = 1; */
/* 	for(i = 1; i <= 2*m-1; i += 2) { /\* (2m-1)!! *\/ */
/*       fact *= i; */
/* 	} */
/* 	pmmint /= fact; */
	
/*   } */
	
/*   /\* pmp1mint = \int_0^\pi P_{m+1}^m(\cos(\theta))\cos(j\theta) d\theta *\/ */


/*   /\*set nonzero elements*\/ */
/*   gsl_matrix_set (m, 0, 1, 0.5); */
/*   gsl_matrix_set (m, 1, 0, 1.0); */
/*   gsl_matrix_set (m, 1, 2, 0.5); */
/*   for(i=2; i<n; i++){ */
/*     gsl_matrix_set (m, i, i-1, 0.5); */
/*     if(i+1<n) */
/*       gsl_matrix_set (m, i, i+1, 0.5); */
/*   }  */
/* } */

/*********************************************************************/
/* Evaluate the integral:                                            */
/*  \int_0^\pi P_m^m(\cos(\theta))\cos(j\theta)\sin(\theta) d\theta  */
/*********************************************************************/
double int_pmm_cosjt_sint(int m, int j)
{
  double pmmint;
  double fact;
  int i;

  if(m%2)
	printf("m is not even in int_pmm_cosjt_sint\n");

  /* j and -j give same result */
  if(j < 0)
	j = -j;

  if(j%2) /* j is odd */
	return 0.0;

  if(j < m) {
	
	pmmint = 2.0*neg1toi(j/2);
	/* (2m-1)!!/(m+j+1)!! */
	fact = 1;
	for(i = 2*m-1; i >= m+j+3; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;

	/* (m+1)!/(m-j+1)!! (next 2 loops) */
	fact = 1;
	for(i = m+1; i >= m-j+2; i--) {
      fact *= i;
	}	
	pmmint *= fact;
	
	fact = 1;
	for(i = m-j; i >= 2; i -= 2) {
      fact *= i;
	}	
	pmmint *= fact;
	
	return pmmint;

  } else if(j == m) {
	
	pmmint = 2.0*neg1toi(j/2); 
	/* (m+1)! */
	fact = 1;
	for(i = m+1; i >= 2; i--) {
      fact *= i;
	}
	pmmint *= fact/(2.0*m+1.0);
	return pmmint;
	
  } else if(j == m+2) {
	
	pmmint = 2.0*neg1toi(j/2); 
	/* (m+1)! */
	fact = 1;
	for(i = m+1; i >= 2; i--) {
      fact *= i;
	}
	pmmint *= fact/((2.0*m+3.0)*(2.0*m+1.0));
	return pmmint;

  }	else { /* j == m+4, m+6, ... */
	
	pmmint = 2.0*neg1toi(j-1-m/2);
	/* divide by (m+j+1)(m+j-1)...(2m+1) */
	fact = 1;
	for(i = m+j+1; i >= 2*m+1; i -= 2) { 
      fact *= i;
	}
	pmmint /= fact;
	/* (m+1)! */
	fact = 1;
	for(i = m+1; i >= 2; i--) { 
      pmmint *= i;
	}
	pmmint *= fact;
	/* (2m-1)!! */
	fact = 1;
	for(i = j-m-3; i >= 2; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;
	return pmmint;
	
  }
}

/***********************************************************************************/
/* Evaluate the integral:                                                          */
/*  \int_0^\pi \cos(\theta)P_{m+1)^m(\cos(\theta))\cos(j\theta)\sin(theta)d\theta  */
/***********************************************************************************/
double int_cost_pmp1m_cosjt_sint(int m, int j)
{
  if(m%2)
	printf("m is not even in int_cost_pmp1m_cosjt_sint\n");
  
  /* j and -j give same result */
  if(j < 0)
	j = -j;
  
  if(j%2) /* j is odd */
	return 0.0;
  else /* j is even */
	return (2.0*m+1.0)*( 0.25*int_pmm_cosjt_sint(m, j+2)
						 + 0.5*int_pmm_cosjt_sint(m, j)
						 + 0.25*int_pmm_cosjt_sint(m, j-2) );
}
