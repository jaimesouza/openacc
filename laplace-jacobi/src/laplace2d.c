#include <stdio.h>
#include <stdlib.h>
#include "sys/time.h"

// return the maximum value
#pragma acc routine vector
int max(float a, float b){
    if( a > b )
        return a;

    return b;
}

int main(){

    // size of grid (width)
    const int N = 4096;

    // size of grid (height)
    const int M = 4096;

    // number of maximum iterations
    const int ITER_MAX = 2000;

    // threshold of convergence
    const float CONV_THRESHOLD = 1.0e-5f;

    // matrix to be solved
    float A[N][M];

    // auxiliary matrix
    float Anew[N][M];

    // seed for random generator
    srand(10);

    // initialize the matrix
    #pragma acc parallel loop
    for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            A[i][j] = rand() % 100;
        }
    }

    float err = 1;
    int iter = 0;

    printf("Jacobi relaxation Calculation: %d x %d mesh\n", N, M);

    struct timeval start_kernel, end_kernel;
    gettimeofday(&start_kernel, NULL);

    /**
     * Jacobi iteration
     * This loop will end if either the maximum change reaches below a set threshold (convergence) or a fixed number of iterations have completed
     */
     #pragma acc data create(Anew[:N][:M]) copy(A[:N][:M])
     while ( err > CONV_THRESHOLD && iter < ITER_MAX ) {

        err = 0.0;

        // calculates the Laplace equation to determine each cell's next value
        #pragma acc parallel loop reduction(max:err) collapse(2)
        for( int i = 1; i < N-1; i++) {
            for(int j = 1; j < M-1; j++) {
                Anew[i][j] = 0.25 * (A[i][j+1] + A[i][j-1] + A[i - 1][j] + A[i+1][j]);
                err = max(err, abs(Anew[i][j] - A[i][j]));
            }
        }

        // copies the next values into the working array for the next iteration
        #pragma acc parallel loop collapse(2)
        for( int i = 1; i < N-1; i++) {
            for( int j = 1; j < M-1; j++ ) {
                A[i][j] = Anew[i][j];
            }
        }

        if(iter % 100 == 0)
            printf("%5d, %0.6f\n", iter, error);

        iter++;
    }

    gettimeofday(&end_kernel, NULL);
    double kernel_time = (double)(end_kernel.tv_sec-start_kernel.tv_sec)+(double)(end_kernel.tv_usec-start_kernel.tv_usec)/1000000;

    printf("\nTime to process the kernel: %f seconds\n", kernel_time );

    return 0;
}
