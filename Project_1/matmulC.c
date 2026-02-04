/* Some stadard libraries for math, data I/O, the time function.. */
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <time.h>
/* Base code for matrix-matrix multiplication. */

double random_matrix(int n,double A[n][n]){
	// This function fills an array with pseudo-random numbers. Operates "in place", no return argument needed.

	int row,column;
	int index=0;
	srand(clock()); // Seed the PRNG (using the number of "clock ticks" since the execution started as seed).

	for(row=0; row<n; row++){
		for(column=0;column<n; column++){
		  A[row][column]=(float)rand() / (float)RAND_MAX; // Random number scaled by its max value.
		};
	}
	  }

double mat_mat(int n, double A[n][n], double B[n][n], double C[n][n]){
  // This is the actual computation of the product:
  int i,j,k;

  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      C[i][j] = 0.0;
      for(k=0; k<n; k++){
	C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}
			   
int main(int argc, char *argv[]){
    int n = atoi(argv[1]); // Get matrix size from command line argument.
    int i,j;
      double (*A)[n] = malloc(sizeof(double[n][n]));
    double (*B)[n] = malloc(sizeof(double[n][n]));
    double (*C)[n] = malloc(sizeof(double[n][n]));
    clock_t time;


    random_matrix(n, A); // Get random matrices.
    random_matrix(n, B);

    time = clock();
    mat_mat(n, A, B, C);
    time = clock() - time;

    printf("%f\n", (double)time/CLOCKS_PER_SEC);


}