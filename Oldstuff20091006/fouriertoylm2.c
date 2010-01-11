/* To compile type: gcc -I/opt/local/include -L/opt/local/lib -lm -lfftw3 -lgsl -lgslcblas fouriertoylm2.c print.c coefficients.c poisson.h */

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
double int_pmm_sinjt_sint(int m, int j);
double int_pmm_cosjt(int m, int j);
double int_pmm_sinjt(int m, int j);

void set_fouriertoylm_meven(gsl_matrix *transform, int m, int nt);
void set_fouriertoylm_modd(gsl_matrix *transform, int m, int nt);
void set_ylmtofourier_meven(gsl_matrix *transform, int m, int nt);
void set_ylmtofourier_modd(gsl_matrix *transform, int m, int nt);

void transform_fouriertoylm(coeff *c_fourier, ylm_coeff *c_ylm, gsl_matrix **fouriertoylm);
void transform_ylmtofourier(ylm_coeff *c_ylm, coeff *c_fourier, gsl_matrix **ylmtofourier);

int main (void)
{
  int z;
  int i;
  int l;
  int m;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  int j;
  gsl_matrix *transform, *transformodd, *inverse, *inverseodd;
  gsl_matrix **fouriertoylm;
  gsl_matrix **ylmtofourier;
  ylm_coeff *c_ylm;
  coeff *c_fourier;

/* for( j = 0; j <= 8; j++) { */
/* 	for(m = 1; m <= 8; m += 2) { */
/* 	  printf("%.18e  ", int_pmm_sinjt(m, j)); */
/* 	} */
/* 	printf("\n"); */
/*   } */
/*   printf("\n"); */
  
/*   m = 4; */
/*   nt = 9; */
/*   transform = gsl_matrix_alloc(nt-m, nt); */
/*   set_fouriertoylm_meven(transform, m, nt); */
/*   print_matrix(transform); */

/*   m = 101; */
/*   nt = 102; */
/*   transformodd = gsl_matrix_alloc(nt-m, nt); */
/*   set_fouriertoylm_modd(transformodd, m, nt);   */
/*   print_matrix(transformodd); */

/*   m = 6; */
/*   nt = 9; */
/*   inverse = gsl_matrix_alloc(nt, nt-m); */
/*   set_ylmtofourier_meven(inverse, m, nt);   */
/*   print_matrix(inverse); */
 
/*   m = 3; */
/*   nt = 9; */
/*   inverseodd = gsl_matrix_alloc(nt, nt-m); */
/*   set_ylmtofourier_modd(inverseodd, m, nt);   */
/*   print_matrix(inverseodd); */

  nz = 1;
  nr = 1;
  nt = 5;
  np = 8; 
  
  npc = np/2 + 1;

  /*>>>>>>>>>>>>>>>> make Fourier to Y_lm matrices <<<<<<<<<<<<<<<*/

  fouriertoylm = malloc(sizeof(gsl_matrix *)*npc);
  /* set even matrices */
  for(m=0; m<nt; m += 2) {
	fouriertoylm[m] = gsl_matrix_alloc(nt-m, nt);
	set_fouriertoylm_meven(fouriertoylm[m], m, nt);
  }
  /* set odd matrices */
  for(m=1; m<nt; m += 2) {
	fouriertoylm[m] = gsl_matrix_alloc(nt-m, nt);
	set_fouriertoylm_modd(fouriertoylm[m], m, nt);
  }

  for(m=0; m<nt; m++) {
	printf("m = %d:\n", m);
	print_matrix(fouriertoylm[m]);
  }

  /*>>>>>>>>>>>>>>>> make Y_lm to Fourier matrices <<<<<<<<<<<<<<<*/

  ylmtofourier = malloc(sizeof(gsl_matrix *)*np);
  /* set even matrices */
  for(m=0; m<nt; m += 2) {
	ylmtofourier[m] = gsl_matrix_alloc(nt, nt-m);
	set_ylmtofourier_meven(ylmtofourier[m], m, nt);
  }
  /* set odd matrices */
  for(m=1; m<nt; m += 2) {
	ylmtofourier[m] = gsl_matrix_alloc(nt, nt-m);
	set_ylmtofourier_modd(ylmtofourier[m], m, nt);
  }

  for(m=0; m<nt; m++) {
	printf("m = %d:\n", m);
	print_matrix(ylmtofourier[m]);
  }

  /*>>>>>>>>>>>>>>>>>> do transforms <<<<<<<<<<<<<<<<<<<<<*/

  c_ylm = ylm_coeff_alloc(nz, nr, nt, npc);
  c_fourier = coeff_alloc(nz, nr, nt, np);

/*   /\* Assign the data. *\/ */
/*   for ( z = 0; z < nz; z++ ) { */
/* 	for ( i = 0; i < nr; i++ ) { */
/* 	  for ( l = 0; l < nt; l++ ) { */
/* 		ylm_coeff_set(c_ylm, z, i, l, 0, REAL, z + i*0.1 + l*0.01); */
/* 		ylm_coeff_set(c_ylm, z, i, l, 0, IMAG, 0.0); */
/* 		for ( m = 1; m <= l; m++ ) { */
/* 		  ylm_coeff_set(c_ylm, z, i, l, m, REAL, z + i*0.1 + l*0.01 + m*0.001); */
/* 		  ylm_coeff_set(c_ylm, z, i, l, m, IMAG, z + i*0.1 + l*0.01 + m*0.001); */
/* 		} */
/* 	  } */
/* 	} */
/*   } */
  
  /* Assign the data. */
  for ( z = 0; z < nz; z++ ) {
	for ( i = 0; i < nr; i++ ) {
	  for ( l = 0; l < nt; l++ ) {
		ylm_coeff_set(c_ylm, z, i, l, 0, REAL, 1.111111);
		ylm_coeff_set(c_ylm, z, i, l, 0, IMAG, 0.0);
		for ( m = 1; m <= l; m++ ) {
		  ylm_coeff_set(c_ylm, z, i, l, m, REAL, 1.111111);
		  ylm_coeff_set(c_ylm, z, i, l, m, IMAG, 1.111111);
		}
		ylm_coeff_set(c_ylm, z, i, nt-1, nt-1, IMAG, 0.0);
	  }
	}
  }
  
  print_ylm_coeff(c_ylm);

  transform_ylmtofourier(c_ylm, c_fourier, ylmtofourier);
  
  print_coeff(c_fourier);
  
  transform_fouriertoylm(c_fourier, c_ylm, fouriertoylm);
  
  print_ylm_coeff(c_ylm);


  /* free matrices as well as the array of pointers to the matrices */
  for(m=0; m<nt; m++) {
	gsl_matrix_free(fouriertoylm[m]);
  }
  free(fouriertoylm);

  /* free matrices as well as the array of pointers to the matrices */
  for(m=0; m<nt; m++) {
	gsl_matrix_free(ylmtofourier[m]);
  }
  free(ylmtofourier);


/*   gsl_matrix_free(transform); */
/*   gsl_matrix_free(transformodd); */
/*   gsl_matrix_free(inverse); */
/*   gsl_matrix_free(inverseodd); */
 
  return 0;
}


void make_fouriertoylm_matrices(gsl_matrix **fouriertoylm, int nt, int npc)
{
  
}


/***************************************************/
/* Transform coefficients in Fourier basis to      */
/* coefficients in spherical harmonic basis.       */
/***************************************************/
void transform_fouriertoylm(coeff *c_fourier, ylm_coeff *c_ylm, gsl_matrix **fouriertoylm)
{
  int z;
  int i;
  int j;
  int m;
  int comp;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  gsl_vector *fourierin;
  gsl_vector *ylmout;

  nz = c_fourier->nz;
  nr = c_fourier->nr;
  nt = c_fourier->nt;
  np = c_fourier->np;

  npc = ( np / 2 ) + 1;

  /* fourierin has same size for all m */
  fourierin = gsl_vector_alloc(nt);
  
  for ( m = 0; m < npc; m++ ) { /* m and k are the same thing */ 
	
	ylmout = gsl_vector_alloc(nt-m);
	
	for(comp=0; comp<2; comp++) {
	  for(z=0; z < nz; z++) {
		for ( i = 0; i < nr; i++ ) {
		  /* fill vector with fourier coefficients */
		  for ( j = 0; j < nt; j++ ) {
			gsl_vector_set(fourierin, j, coeff_get(c_fourier, z, i, j, m, comp));
		  }
		  /* convert to spherical harmonic coefficients */
		  /* using matrix transformation */
		  gsl_blas_dgemv(CblasNoTrans,
						 1.0, fouriertoylm[m], fourierin,
						 0.0, ylmout);
		  /* put spherical harmonic coefficients in new structure */
		  /* Y_lm where j=l<m is zero */
		  for ( j = m; j < nt; j++ ) {
			ylm_coeff_set(c_ylm, z, i, j, m, comp, gsl_vector_get(ylmout, j-m));
		  }
		}
	  }
	}
	gsl_vector_free(ylmout);
  }

  gsl_vector_free(fourierin);

}


/**********************************************************/
/* Transform coefficients in spherical harmonic basis to  */
/* coefficients in Fourier basis.                         */
/**********************************************************/
void transform_ylmtofourier(ylm_coeff *c_ylm, coeff *c_fourier, gsl_matrix **ylmtofourier)
{
  int z;
  int i;
  int j;
  int m;
  int comp;
  int nz;
  int nr;
  int nt;
  int np;
  int npc;
  gsl_vector *ylmin;
  gsl_vector *fourierout;
  
  nz = c_fourier->nz;
  nr = c_fourier->nr;
  nt = c_fourier->nt;
  np = c_fourier->np;
  
  npc = ( np / 2 ) + 1;

  /* fourierout has same size for all m */
  fourierout = gsl_vector_alloc(nt);
  
  for ( m = 0; m < npc; m++ ) { /* m and k are the same thing */ 
	/*printf("m=%d:\n", m);*/
	ylmin = gsl_vector_alloc(nt-m);
	
	for(comp=0; comp<2; comp++) {
	  for(z=0; z < nz; z++) {
		for ( i = 0; i < nr; i++ ) {
		  /* fill vector with spherical harmonic coefficients */
		  /* Y_lm where j=l<m is zero */
		  for ( j = m; j < nt; j++ ) {
			/*printf("%d, %d, %d, %d, %d, %f\n", z, i, j, m, comp, 
			  ylm_coeff_get(c_ylm, z, i, j, m, comp));*/
			gsl_vector_set(ylmin, j-m, ylm_coeff_get(c_ylm, z, i, j, m, comp));
			/*print_vector(ylmin);*/
		  }
		  /* convert to fourier coefficients */
		  /* using matrix transformation */
		  gsl_blas_dgemv(CblasNoTrans,
						 1.0, ylmtofourier[m], ylmin,
						 0.0, fourierout);
		  /*print_vector(fourierout);*/
		  /* put Fourier coefficients in new structure */
		  for ( j = 0; j < nt; j++ ) {
			coeff_set(c_fourier, z, i, j, m, comp, gsl_vector_get(fourierout, j));
		  }
		}
	  }
	}
	gsl_vector_free(ylmin);
  }

  gsl_vector_free(fourierout);

}



/*****************************************************************/
/* Make matrix for converting cos(j theta) to P_l^m(cos(theta))  */
/* where m is even.                                              */
/*****************************************************************/
void set_fouriertoylm_meven(gsl_matrix *transform, int m, int nt)
{
  int row; /* number of rows / length of output vector */
  int col; /* number of columns */
  int i, j;
  int l;
  gsl_matrix *Imatrix;
  int n; /* for calculating factorials */
  double fact; 
  double normalized;
  
  row = transform->size1; 
  col = transform->size2; 

  /* make sure dimensions are correct */
  if(m%2) {
	printf("m is not even in set_fouriertoylm_meven\n");
  } else if(nt != col) {
	printf("nt does not match number of columns in transform matrix ");
	printf("in set_fouriertoylm_meven\n");
  } else if(m != nt - row) {
	printf("nt-m does not match number of rows in transform matrix ");
	printf("in set_fouriertoylm_meven\n");
  } else if(m >= nt) {
	printf("m cannot be >= nt in set_fouriertoylm_meven\n");
  }

  /* allocate memory and initialize intermediate matrices */
  Imatrix = gsl_matrix_alloc(nt, 2*nt-m-1);
  gsl_matrix_set_zero(Imatrix);

    
  /* fill mth row of I (only even j are non zero) */
  for( j = 0; j < 2*nt-m-1; j += 2) {
	gsl_matrix_set(Imatrix, m, j, int_pmm_cosjt_sint(m, j));
  }
  
  if (row > 1) {
	/* fill m+1 row of I (only odd j are non zero) */
	for( j = 1; j < 2*nt-m-1; j += 2) {
	  gsl_matrix_set(Imatrix, m+1, j,
					 0.5*(2.0*m+1.0)*(int_pmm_cosjt_sint(m, j-1) + int_pmm_cosjt_sint(m, j+1)));
	}
  }
  
  /* fill the other rows of I */
  /* reflect indices of matrix if needed */
  /* ex: j = 0-2 goes to j = 2 */
  for( l = m+2; l < nt; l++) {
	for( j = 0; j < 2*nt-m-3; j++) {
	  gsl_matrix_set(Imatrix, l, j, 
					 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix, l-1, ABS(j-1)) 
					 + 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix, l-1, j+1)
					 - (l+m-1.0)/(l-m)*gsl_matrix_get(Imatrix, l-2, j));	 
	}
  }
  
  /*print_matrix(Imatrix);*/
  
  
  /*initialize matrix*/
  gsl_matrix_set_zero(transform);
 
  /*>>>>>>>>>>> multiply by normalization and plm to ylm factors <<<<<<<<<*/

  /* Fill even rows */
  for( i = 0; i < row; i += 2 ) {
	for( j = 0; j < col; j += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = sqrt(PI*(2*l+1)/fact)*gsl_matrix_get(Imatrix, l, j);
	  gsl_matrix_set(transform, i, j, normalized);
	}
  }
  /* Fill odd rows */
  for( i = 1; i < row; i += 2 ) {
	for( j = 1; j < col; j += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = sqrt(PI*(2*l+1)/fact)*gsl_matrix_get(Imatrix, l, j);
	  gsl_matrix_set(transform, i, j, normalized);
	}
  }

  gsl_matrix_free(Imatrix);
}


/*****************************************************************/
/* Make matrix for converting sin(j theta) to P_l^m(cos(theta))  */
/* where m is odd.                                               */
/*****************************************************************/
void set_fouriertoylm_modd(gsl_matrix *transform, int m, int nt)
{
  int row; /* number of rows / length of output vector */
  int col; /* number of columns */
  int i, j;
  int l;
  gsl_matrix *Imatrix; 
  int n; /* for calculating factorials */
  double fact; 
  double normalized;
  
  row = transform->size1; 
  col = transform->size2; 

  /* make sure dimensions are correct */
  if(!(m%2)) {
	printf("m is not odd in set_fouriertoylm_modd\n");
  } else if(nt != col) {
	printf("nt does not match number of columns in transform matrix ");
	printf("in set_fouriertoylm_modd\n");
  } else if(m != nt - row) {
	printf("nt-m does not match number of rows in transform matrix ");
	printf("in set_fouriertoylm_modd\n");
  } else if(m >= nt) {
	printf("m cannot be >= nt in set_fouriertoylm_modd\n");
  }

  /* allocate memory and initialize intermediate matrix to 0 */
  Imatrix = gsl_matrix_alloc(nt, 2*nt-m-1);
  gsl_matrix_set_zero(Imatrix);

    
  /* fill mth row of I (only odd j are non zero) */
  for( j = 1; j < 2*nt-m-1; j += 2) {
	gsl_matrix_set(Imatrix, m, j, int_pmm_sinjt_sint(m, j));
  }
  
  if (row > 1) {
	/* fill m+1 row of I (only even j are non zero (also 0 for j=0)) */
	for( j = 2; j < 2*nt-m-1; j += 2) {
	  gsl_matrix_set(Imatrix, m+1, j,
					 0.5*(2.0*m+1.0)*(int_pmm_sinjt_sint(m, j-1) + int_pmm_sinjt_sint(m, j+1)));
	}
  }
  
  /* fill the other rows of I */
  /* start at j=1 because all j=0 elements are zero */
  for( l = m+2; l < nt; l++) {
	for( j = 1; j < 2*nt-m-3; j++) {
	  gsl_matrix_set(Imatrix, l, j, 
					 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix, l-1, j-1) 
					 + 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix, l-1, j+1)
					 - (l+m-1.0)/(l-m)*gsl_matrix_get(Imatrix, l-2, j));	 
	}
  }
  
  /*print_matrix(Imatrix);*/
  
  
  /*initialize matrix*/
  gsl_matrix_set_zero(transform);
 
  /*>>>>>>>>>>> multiply by normalization and plm to ylm factors <<<<<<<<<*/

  /* Fill even rows (odd elements are nonzero) */
  for( i = 0; i < row; i += 2 ) {
	for( j = 1; j < col; j += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = sqrt(PI*(2*l+1)/fact)*gsl_matrix_get(Imatrix, l, j);
	  gsl_matrix_set(transform, i, j, normalized);
	}
  }
  /* Fill odd rows (even nonzero elements are nonzero */
  for( i = 1; i < row; i += 2 ) {
	for( j = 2; j < col; j += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = sqrt(PI*(2*l+1)/fact)*gsl_matrix_get(Imatrix, l, j);
	  gsl_matrix_set(transform, i, j, normalized);
	}
  }

  gsl_matrix_free(Imatrix);
}


/*****************************************************************/
/* Make matrix for converting P_l^m(cos(theta)) to cos(j theta)  */
/* where m is even.                                              */
/*****************************************************************/
void set_ylmtofourier_meven(gsl_matrix *transform, int m, int nt)
{
  int row; /* number of rows / length of output vector */
  int col; /* number of columns */
  int i, j;
  int l;
  gsl_matrix *Imatrix; 
  int n; /* for calculating factorials */
  double fact; 
  double normalized;

  row = transform->size1; 
  col = transform->size2; 

  /* make sure dimensions are correct */
  if(m%2) {
	printf("m is not even in set_ylmtofourier_meven\n");
  } else if(nt != row) {
	printf("nt does not match number of rows in transform matrix ");
	printf("in set_ylmtofourier_meven\n");
  } else if(m != nt - col) {
	printf("nt-m does not match number of columns in transform matrix ");
	printf("in set_ylmtofourier_meven\n");
  } else if(m >= nt) {
	printf("m cannot be >= nt in set_ylmtofourier_meven\n");
  }

  /* allocate memory and initialize intermediate matrices */
  Imatrix = gsl_matrix_alloc(2*nt-m-1, nt);
  gsl_matrix_set_zero(Imatrix);

    
  /* fill mth column of I (only even j are non zero) */
  for( j = 0; j < 2*nt-m-1; j += 2) {
	gsl_matrix_set(Imatrix, j, m, int_pmm_cosjt(m, j));
  }
  
  if (col > 1) {
	/* fill m+1 column of I (only odd j are non zero) */
	for( j = 1; j < 2*nt-m-1; j += 2) {
	  gsl_matrix_set(Imatrix, j, m+1,
					 0.5*(2.0*m+1.0)*(int_pmm_cosjt(m, j-1) + int_pmm_cosjt(m, j+1)));
	}
  }
  
  /*print_matrix(Imatrix);*/

  /* fill the other columns of I */
  /* reflect indices of matrix if needed */
  /* ex: j = 0-1 goes to j = +1 */
  for( l = m+2; l < nt; l++) {
	for( j = 0; j < 2*nt-m-3; j++) {
	  gsl_matrix_set(Imatrix, j, l, 
					 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix, ABS(j-1), l-1) 
					 + 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix,  j+1, l-1)
					 - (l+m-1.0)/(l-m)*gsl_matrix_get(Imatrix, j, l-2));
	  /*printf("j=%d, l=%d, %f\n", j, l, gsl_matrix_get(Imatrix, j, l));*/
	}
	/*printf("\n");*/
  }
  
  /*print_matrix(Imatrix);*/
  
  
  /*initialize matrix*/
  gsl_matrix_set_zero(transform);
 
  /*>>>>>>>>>>> multiply by normalization and plm to ylm factors <<<<<<<<<*/

  /* Fill even columns */
  for( j = 0; j < row; j += 2 ) {
	for( i = 0; i < col; i += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = (2-delta(0, j))/PI*sqrt((2*l+1)/(4*PI*fact))*gsl_matrix_get(Imatrix, j, l);
	  gsl_matrix_set(transform, j, i, normalized);
	}
  }
  /* Fill odd columns */
  for( j = 1; j < row; j += 2 ) {
	for( i = 1; i < col; i += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = (2-delta(0, j))/PI*sqrt((2*l+1)/(4*PI*fact))*gsl_matrix_get(Imatrix, j, l);
	  gsl_matrix_set(transform, j, i, normalized);
	}
  }

  gsl_matrix_free(Imatrix);
}



/*****************************************************************/
/* Make matrix for converting P_l^m(cos(theta)) to sin(j theta)  */
/* where m is odd.                                               */
/*****************************************************************/
void set_ylmtofourier_modd(gsl_matrix *transform, int m, int nt)
{
  int row; /* number of rows / length of output vector */
  int col; /* number of columns */
  int i, j;
  int l;
  gsl_matrix *Imatrix; 
  int n; /* for calculating factorials */
  double fact; 
  double normalized;

  row = transform->size1; 
  col = transform->size2; 

  /* make sure dimensions are correct */
  if(!(m%2)) {
	printf("m is not odd in set_ylmtofourier_modd\n");
  } else if(nt != row) {
	printf("nt does not match number of rows in transform matrix ");
	printf("in set_ylmtofourier_modd\n");
  } else if(m != nt - col) {
	printf("nt-m does not match number of columns in transform matrix ");
	printf("in set_ylmtofourier_modd\n");
  } else if(m >= nt) {
	printf("m cannot be >= nt in set_ylmtofourier_modd\n");
  }

  /* allocate memory and initialize intermediate matrices */
  Imatrix = gsl_matrix_alloc(2*nt-m-1, nt);
  gsl_matrix_set_zero(Imatrix);

    
  /* fill mth column of I (only odd j are nonzero) */
  for( j = 1; j < 2*nt-m-1; j += 2) {
	gsl_matrix_set(Imatrix, j, m, int_pmm_sinjt(m, j));
  }
  
  if (col > 1) {
	/* fill m+1 column of I (only even, nonzero j are nonzero) */
	for( j = 2; j < 2*nt-m-1; j += 2) {
	  gsl_matrix_set(Imatrix, j, m+1,
					 0.5*(2.0*m+1.0)*(int_pmm_sinjt(m, j-1) + int_pmm_sinjt(m, j+1)));
	}
  }
  
  /*print_matrix(Imatrix);*/

  /* fill the other columns of I */
  /* first row is zero so start j at 1 */
  for( l = m+2; l < nt; l++) {
	for( j = 1; j < 2*nt-m-3; j++) {
	  gsl_matrix_set(Imatrix, j, l, 
					 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix, j-1, l-1) 
					 + 0.5*(2.0*l-1.0)/(l-m)*gsl_matrix_get(Imatrix,  j+1, l-1)
					 - (l+m-1.0)/(l-m)*gsl_matrix_get(Imatrix, j, l-2));
	  /*printf("j=%d, l=%d, %f\n", j, l, gsl_matrix_get(Imatrix, j, l));*/
	}
	/*printf("\n");*/
  }
  
  /*print_matrix(Imatrix);*/
  
  
  /*initialize matrix*/
  gsl_matrix_set_zero(transform);
 
  /*>>>>>>>>>>> multiply by normalization and plm to ylm factors <<<<<<<<<*/

  /* Fill even columns */
  for( j = 0; j < row; j += 2 ) {
	for( i = 1; i < col; i += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = (2-delta(0, j))/PI*sqrt((2*l+1)/(4*PI*fact))*gsl_matrix_get(Imatrix, j, l);
	  gsl_matrix_set(transform, j, i, normalized);
	}
  }
  /* Fill odd columns */
  for( j = 1; j < row; j += 2 ) {
	for( i = 0; i < col; i += 2 ) {
	  l=i+m;
	  /* (l+m)!\(l-m)! */
	  fact = 1;
	  for( n = l+m; n >= l-m+1; n-- ) {
		fact *= n;
	  }
	  normalized = (2-delta(0, j))/PI*sqrt((2*l+1)/(4*PI*fact))*gsl_matrix_get(Imatrix, j, l);
	  gsl_matrix_set(transform, j, i, normalized);
	}
  }
  
  gsl_matrix_free(Imatrix);
}





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


/*********************************************************************/
/* Evaluate the integral:                                            */
/*  \int_0^\pi P_m^m(\cos(\theta))\sin(j\theta)\sin(\theta) d\theta  */
/*********************************************************************/
double int_pmm_sinjt_sint(int m, int j)
{
  double pmmint;
  double fact;
  int i;
  int jsign = 1; /* later switched to -1 if j is negative */
  
  if(!(m%2))
	printf("m is not odd in int_pmm_sinjt_sint\n");
  
  /* j and -j give results that differ by a sign */
  if(j < 0){
	j = -j; /* switch the sign later in the return statement*/
	jsign = -1;
  }
  
  if(!(j%2)) /* j is even */
	return 0.0;
  
  if(j < m) {
	
	pmmint = neg1toi((j+1)/2)*pow(2.0, (m+3)/2);
	/* ((m+1)/2)! */
	fact = 1;
	for(i = (m+1)/2; i >= 2; i--) {
      fact *= i;
	}
	pmmint *= fact;

	/* (2m-1)!!/(m+j+1)!! */
	fact = 1;
	for(i = 2*m-1; i >= m+j+3; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;
	
	/* m!!/(m-j+1)!! */
	fact = 1;
	for(i = m; i >= m-j+3; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;
	
	return jsign*pmmint;

  } else if(j == m) {
	
	pmmint = neg1toi((j+1)/2)*pow(2.0, (m+3)/2);
	/* ((m+1)/2)! */
	fact = 1;
	for(i = (m+1)/2; i >= 2; i--) {
      fact *= i;
	}
	pmmint *= fact;
	
	/* m!! */
	fact = 1;
	for(i = m; i >= 2; i -= 2) {
      fact *= i;
	}
	pmmint *= fact/(2.0*m+1.0);

	return jsign*pmmint;
	
  } else if(j == m+2) {
	
	pmmint = neg1toi((j+1)/2)*pow(2.0, (m+3)/2);
	/* ((m+1)/2)! */
	fact = 1;
	for(i = (m+1)/2; i >= 2; i--) {
      fact *= i;
	}
	pmmint *= fact;
	
	/* m!! */
	fact = 1;
	for(i = m; i >= 2; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;
	
	return jsign*pmmint/((2.0*m+3.0)*(2.0*m+1.0));

  }	else { /* j == m+4, m+6, ... */
	
	pmmint = neg1toi((m+3)/2)*pow(2.0, (m+3)/2);

	/* ((m+1)/2)! */
	fact = 1;
	for(i = (m+1)/2; i >= 2; i--) {
      fact *= i;
	}
	pmmint *= fact;
	
	/* m!! */
	fact = 1;
	for(i = m; i >= 2; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;
	
	/* (j-m-3)!! */
	fact = 1;
	for(i = j-m-3; i >= 2; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;

	/* divide by (m+j+1)(m+j-1)...(2m+1) */
	fact = 1;
	for(i = m+j+1; i >= 2*m+1; i -= 2) { 
      fact *= i;
	}
	pmmint /= fact;

	return jsign*pmmint;
	
  }
}


/*********************************************************************/
/* Evaluate the integral:                                            */
/*  \int_0^\pi P_m^m(\cos(\theta))\cos(j\theta) d\theta              */
/*********************************************************************/
double int_pmm_cosjt(int m, int j)
{
  double pmmint;
  double fact;
  int i;
  
  if(m%2)
	printf("m is not even in int_pmm_cosjt\n");
  
  /* j and -j give same results */
  if(j < 0)
	j = -j;
  
  if(j%2) /* j is odd */
	return 0.0;
  else if(j > m)
	return 0.0;
  else {
	pmmint = neg1toi(j/2) * PI / pow(2.0, m);

	/* multiply by m! / ((m+j)/2)! */
	fact = 1;
	for(i = m; i >= (m+j+2)/2; i--) {
      fact *= i;
	}
	pmmint *= fact;
	
	/* multiply by (2m-1)!! */
	fact = 1;
	for(i = 2*m-1; i >= 2; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;	

	/* divide by ((m-j)/2)! */
	fact = 1;
	for(i = (m-j)/2; i >= 2; i--) {
      fact *= i;
	}
	pmmint /= fact;

	return pmmint;
  }

}


/*********************************************************************/
/* Evaluate the integral:                                            */
/*  \int_0^\pi P_m^m(\cos(\theta))\sin(j\theta) d\theta              */
/*********************************************************************/
double int_pmm_sinjt(int m, int j)
{
  double pmmint;
  double fact;
  int i;
  int jsign = 1; /* later switched to -1 if j is negative */

  if(!(m%2))
	printf("m is not odd in int_pmm_sinjt\n");
  
  /* j and -j give results that differ by a sign */
  if(j < 0){
	j = -j; /* switch the sign later in the return statement*/
	jsign = -1;
  }

  if(!(j%2)) /* j is even */
	return 0.0;
  else if(j > m)
	return 0.0;
  else {
	pmmint = neg1toi((j+1)/2) * PI / pow(2.0, m);
	
	/* multiply by m! / ((m+j)/2)! */
	fact = 1;
	for(i = m; i >= (m+j+2)/2; i--) {
      fact *= i;
	}
	pmmint *= fact;
	
	/* multiply by (2m-1)!! */
	fact = 1;
	for(i = 2*m-1; i >= 2; i -= 2) {
      fact *= i;
	}
	pmmint *= fact;	

	/* divide by ((m-j)/2)! */
	fact = 1;
	for(i = (m-j)/2; i >= 2; i--) {
      fact *= i;
	}
	pmmint /= fact;

	return jsign*pmmint;
  }

}






/***********************************************************************************/
/* Evaluate the integral:                                                          */
/*  \int_0^\pi \cos(\theta)P_{m+1)^m(\cos(\theta))\cos(j\theta)\sin(theta)d\theta  */
/***********************************************************************************/
/* double int_cost_pmp1m_cosjt_sint(int m, int j) */
/* { */
/*   if(m%2) */
/* 	printf("m is not even in int_cost_pmp1m_cosjt_sint\n"); */
  
/*   /\* j and -j give same result *\/ */
/*   if(j < 0) */
/* 	j = -j; */
  
/*   if(j%2) /\* j is odd *\/ */
/* 	return 0.0; */
/*   else /\* j is even *\/ */
/* 	return (2.0*m+1.0)*( 0.25*int_pmm_cosjt_sint(m, j+2) */
/* 						 + 0.5*int_pmm_cosjt_sint(m, j) */
/* 						 + 0.25*int_pmm_cosjt_sint(m, j-2) ); */
/* } */

