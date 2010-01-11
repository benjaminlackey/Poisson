#include <stdio.h>
#include <stdlib.h>
void makematrix(float **d, int n);
float **fmatrix( int xlo, int ylo, int xhi, int yhi );
void freefmatrix( float **m, int xlo, int ylo );

int main()
{
  
  int i, j, n;
  float **pointer;
  
  n=4;
  
  pointer = fmatrix(0, 10, 0, 10);

  pointer[2][2] = 7;
  printf("%f\t", pointer[2][2]);

 /*  makematrix(pointer, n); */
  
/*   for(i=0; i<=n; i++){ */
/*     for(j=0; j<=n; j++){ */
/*       printf("%f\t", pointer[i][j]); */
/*     } */
/*     printf("\n"); */
/*   } */
  
  return 0;
}

void makematrix(float **d, int n)
{
  /*initialize matrix*/
  int i, j;
  for(i=0; i<=n; i++){
    for(j=0; j<=n; j++){
      d[i][j] = 3.0;
    }
  }
}


float **fmatrix( int xlo,
                 int ylo,
                 int xhi,
                 int yhi )
/*------------------------------------------------------------------*/
/* PURPOSE: Allocate space for a matrix with subscript range        */
/*          M[ylo..yhi][xlo..xhi] for floats.                       */
/*------------------------------------------------------------------*/
{
   int      rows = yhi - ylo + 1;
   int      cols = xhi - xlo + 1;
   int      i;
   float    **m;


   /* Allocate memory for pointers to rows */
   m = (float **) malloc( rows * sizeof(float*) );
   if (!m)
      return( NULL );
   m -= ylo;                                /* Adjust subscript */

   /* Pointer to first row allocates memory for box */
   m[ylo] = (float *) malloc( rows * cols * sizeof(float) );
   if (!m[ylo])
       return( NULL );
   m[ylo] -= xlo;                           /* Adjust subscript */

   /* Set pointers to rows */
   for (i = ylo+1; i <= yhi; i++)
      m[i] = m[i-1] + cols;

   /* Return pointer to array of pointers to rows */
   return( m );
}



void freefmatrix( float **m,
                  int   xlo,
                  int   ylo )
/*------------------------------------------------------------------*/
/* PURPOSE: Free space allocated by a matrix made with function     */
/*          'fmatrix'.                                              */
/*------------------------------------------------------------------*/
{
   free( m[ylo] + xlo );
   free( m + ylo );
}
