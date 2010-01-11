#include <stdio.h>

#define swap(a, b) tempr=(a); (a)=(b); (b)=tempr

void fft(double data[], unsigned long nn, int isign);

main()
{
  int i;
  int isign = 1;
  unsigned long nn = 8;
  double data[16]={0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
  
  for(i=0; i<2*nn; i++)
	printf("%f\n", data[i]);
  
  fft(data, nn, isign);
  printf("hello\n");
  for(i=0; i<2*nn; i++) 
 	printf("%f\n", data[i]); 
}


/* nn must be integer power of 2 */
void fft(double data[], unsigned long nn, int isign)
{
  /* nn number of elements in array data */
  
  unsigned long n;
  unsigned long i, j; 
  unsigned long mmax; 
  unsigned long m; 
  unsigned long istep;
  double wtemp, wr, wpr, wpi, wi, theta;
  double tempr; /* temporary variable in swap macro */ 
  double tempi;

  n = nn << 1; /* left bit shift (multiply by 2) */
  j=0;
  for(i=0; i<n; i+=2){
	if(j>i){
	  swap(data[j], data[i]);
	  swap(data[j+1], data[i+1]);
	}
	m = n >> 1;
	while(m>=2 && j>m){
	  j -= m;
	  m >>= 1;
	}
	j += m;
  }
}
