/* gcc -g -pedantic -Wall coeff.c */

/* make error checking optional */
/* make functions inline */

#include <stdio.h>
#include <stdlib.h>

typedef struct 
{
  int nz;
  int nn;
  int nl;
  double *data;
} coeff;

typedef struct
{
  int nz;
  int nl;
  double *data;
} bound_coeff;

#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */

coeff *coeff_alloc(int z, int n, int L);
void coeff_set(coeff *c, int z, int n, int L, int m, int i, double x);
double coeff_get(coeff *c, int z, int n, int L, int m, int i);
void print_coeff(coeff *c);

bound_coeff *bound_coeff_alloc(int z, int L);
void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x);
double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i);
void print_bound_coeff(bound_coeff *b);


int main(void)
{
  int z, n, L, m, k;
  int nz = 1;
  int nn = 1;
  int nl = 4;
  coeff *c_znlm = coeff_alloc(nz, nn, nl);
  bound_coeff *b_zlm = bound_coeff_alloc(nz, nl);

  k=0;
  for(z=0; z<nz; z++) 
    for(n=0; n<nn; n++) 
      for(L=0; L<nl; L++){ 
		printf("k=%d\n", k);	
		coeff_set(c_znlm, z, n, L, 0, 0, k); 
		printf("%d\t%d\t%d\t%d\t%d\t%f\n", z, n, L, 0, 0, coeff_get(c_znlm, z, n, L, 0, 0));
		k++;
		for(m=1; m<=L; m++) {
		  printf("k=%d\n", k);
		  coeff_set(c_znlm, z, n, L, m, 0, k); 
		  printf("%d\t%d\t%d\t%d\t%d\t%f\n", z, n, L, m, 0, coeff_get(c_znlm, z, n, L, m, 0));
		  k++;
		  
		  printf("k=%d\n", k);
		  coeff_set(c_znlm, z, n, L, m, 1, k); 
		  printf("%d\t%d\t%d\t%d\t%d\t%f\n", z, n, L, m, 1, coeff_get(c_znlm, z, n, L, m, 1));
		  k++;
		}
	  }
  
  print_coeff(c_znlm);

  
  k=0;
  for(z=0; z<nz; z++)  
	for(L=0; L<nl; L++){ 
	  printf("k=%d\n", k);	
	  bound_coeff_set(b_zlm, z, L, 0, 0, k); 
	  printf("%d\t%d\t%d\t%d\t%f\n", z, L, 0, 0, bound_coeff_get(b_zlm, z, L, 0, 0));
	  k++;
	  for(m=1; m<=L; m++) {
		printf("k=%d\n", k);
		bound_coeff_set(b_zlm, z, L, m, 0, k); 
		printf("%d\t%d\t%d\t%d\t%f\n", z, L, m, 0, bound_coeff_get(b_zlm, z, L, m, 0));
		k++;
		
		printf("k=%d\n", k);
		bound_coeff_set(b_zlm, z, L, m, 1, k); 
		printf("%d\t%d\t%d\t%d\t%f\n", z, L, m, 1, bound_coeff_get(b_zlm, z, L, m, 1));
		k++;
	  }
	}
  
  print_bound_coeff(b_zlm);
  return 0;
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
  /* m goes from -l to l for l */
  
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


/********************************/
/* Set coefficient c_znlm to x. */
/********************************/
/*extern inline */
void coeff_set(coeff *c, int z, int n, int L, int m, int i, double x)
{ 
  /* do error checking */
  if((z >= c->nz)||(n >= c->nn)||(L >= c->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in coeff_set\n");
  
  /* index = z*nn*nl^2 + n*nl^2 + L^2 + 2m + i + (delta_{0m} - 1) */
  /*       = ((z*nn + n)*nl^2 + L^2 + 2m + i + (delta_{0m} - 1) */
  c->data[(z*c->nn + n)*c->nl*c->nl + L*L + 2*m + i + delta(0, m) - 1] = x ;
}


/********************************************/
/* Set coefficient b_zlm for boundary to x. */
/********************************************/
/*extern inline */
void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x)
{ 
  /* do error checking */
  if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_set\n");
  
  /* index = z*nl^2 + n*nl^2 + L^2 + 2m + i+ (delta_0m - 1) */
  b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1] = x ;
}


/*****************************/
/* Get value of c_znlm.      */
/*****************************/
/*extern inline */
double coeff_get(coeff *c, int z, int n, int L, int m, int i)
{
  /* do error checking */
  if((z >= c->nz)||(n >= c->nn)||(L >= c->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in coeff_get\n");
  
  return c->data[(z*c->nn + n)*c->nl*c->nl + L*L + 2*m + i + delta(0, m) - 1];
} 

/*****************************/
/* Get value of b_zlm.      */
/*****************************/
/*extern inline */
double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i)
{
  /* do error checking */
  if((z >= b->nz)||(L >= b->nl)||(m > L)||(m < 0)||(i<0)||(i>1))
    printf("index out of bounds in bound_coeff_get\n");
  
  return b->data[z*b->nl*b->nl + L*L + 2*m + i + delta(0, m) - 1];
}


/*******************************/
/* Print all the coefficients. */
/*******************************/
void print_coeff(coeff *c){
  int nz = c->nz;
  int nn = c->nn;
  int nl = c->nl;
  int z, n, L, m;
  
  for(z=0; z<nz; z++) {
    printf("zone %d:\n", z);
    for(n=0; n<nn; n++) {
      printf("radial basis %d:\n\t", n);
	  printf("m=0\t");
      for(m=1; m<nl; m++)
		printf("m=%d\t\t", m);
      printf("\n");
	  printf("L=0:\t");
	  printf("(%f, -)\n", coeff_get(c, z, n, 0, 0, 0));
      for(L=1; L<nl; L++) {
		printf("L=%d:\t", L);
		printf("(%f, -)\t", coeff_get(c, z, n, L, 0, 0));
		for(m=1; m<=L; m++){
		  printf("(%f, %f)\t", coeff_get(c, z, n, L, m, 0), coeff_get(c, z, n, L, m, 1));
		}
		printf("\n");
      }
    }
  }
  printf("\n");
}


/****************************************/
/* Print all the boundary coefficients. */
/****************************************/
void print_bound_coeff(bound_coeff *b){
  int nz = b->nz;
  int nl = b->nl;
  int z, L, m;
  
  for(z=0; z<nz; z++) {
    printf("zone %d:\n", z);
	printf("\tm=0\t");
	for(m=1; m<nl; m++)
	  printf("m=%d\t\t", m);
	printf("\n");
	printf("L=0:\t");
	printf("(%f, -)\n", bound_coeff_get(b, z, 0, 0, 0));
	for(L=1; L<nl; L++) {
	  printf("L=%d:\t", L);
	  printf("(%f, -)\t", bound_coeff_get(b, z, L, 0, 0));
	  for(m=1; m<=L; m++){
		printf("(%f, %f)\t", bound_coeff_get(b, z, L, m, 0), bound_coeff_get(b, z, L, m, 1));
	  }
	  printf("\n");
	}
  }
  printf("\n");
}

