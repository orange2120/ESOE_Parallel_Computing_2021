#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#define N 8

enum MessageTag {
    EMPTY_MSG,   // The message is for request
    RESULT_MSG,  // The computation result
    N_MSG,       // The exponent n, receive -1 for not perfect number
};

void manager(int rank, int size);
void worker(int rank, int size, MPI_Comm commWorker);

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2) {
        fprintf(stderr, "This program needs at least 2 processes to run\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm commWorker;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &commWorker);
        manager(rank, size);
    } else {
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &commWorker);
        worker(rank, size, commWorker);
    }

    MPI_Finalize();
    return 0;
}

void manager(int rank, int size) {
    int exponent = 2;
    int perfectCnt = 0;
    int nRunning = size - 1;
    int n;

    MPI_Status status;

    double elapsedTime = -MPI_Wtime();
    while (nRunning) {
        MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == RESULT_MSG) {
            if (n > 0) {
                // this is a perfect number
                // TODO
                printf("Perfect Number: %llu\n", ((1ull << n) - 1) * (1 << (n - 1)));
                ++perfectCnt;
            }
        }
        if (perfectCnt < N && exponent < 50) {
            // tell worker to calculate the next exponent
            // printf("Manager pass %d to Worker %d\n", exponent, status.MPI_SOURCE);
            MPI_Send(&exponent, 1, MPI_INT, status.MPI_SOURCE, N_MSG, MPI_COMM_WORLD);
        } else {
            // tell worker to terminate
            MPI_Send(NULL, 0, MPI_INT, status.MPI_SOURCE, N_MSG,
                     MPI_COMM_WORLD);
            --nRunning;
        }

        ++exponent;
    }
    elapsedTime += MPI_Wtime();
    printf("Time elapsed: %.8lf\n", elapsedTime);
    fflush(stdout);
}

int checkPrime(unsigned long long n) {
    for (int i = 3; i <= (int)sqrt(n); i += 2) {
        if (n % i == 0) return 0;
    }
    return 1;
}

void worker(int rank, int size, MPI_Comm commWorker) {
    MPI_Status status;
    int n = 0;
    MPI_Send(NULL, 0, MPI_INT, EMPTY_MSG, 0, MPI_COMM_WORLD);
    while (1) {
        MPI_Recv(&n, 1, MPI_INT, 0, N_MSG, MPI_COMM_WORLD, &status);
        int recvCnt;
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recvCnt);
        if (recvCnt == 0) {
            // printf("Worker %d terminated!\n", rank);
            break;
        }
        unsigned long long num = (1ull << n) - 1;
        // check prime
        if (!checkPrime(num)) n = -1;  // not prime, send -1 back to manager
        MPI_Send(&n, 1, MPI_INT, 0, RESULT_MSG, MPI_COMM_WORLD);
    }
}
