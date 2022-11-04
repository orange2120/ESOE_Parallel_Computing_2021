#define CNT 8
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

/**
 * 1
 * 6
 * 28
 * 496
 * 8128
 * 33550336
 * 8589869056
 * 137438691328
 * 2305843008139952128
 */

int checkPrime(int n) {
    for (int i = 2; i <= sqrt(n); ++i) {
        if (n % i == 0) return 0;
    }
    return 1;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int cnt = 0;
    int n = 1;
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsedTime = -MPI_Wtime();
    while (1) {
        if (checkPrime(((1ll << n) - 1))) {
            printf("%d, %lld\n", n, ((1ll << n) - 1) * (1 << (n - 1)));
            ++cnt;
        }
        ++n;
        if (cnt > CNT) break;
    }
    elapsedTime += MPI_Wtime();
    printf("Time elapsed: %.8lf\n", elapsedTime);
    fflush(stdout);
    MPI_Finalize();
    return 0;
}
