/*
 * Sum of 1+2+...+p 
 */

#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    int id; /* Process rank */
    int p;  /* Number of processes */
    int n;
    int sum_par;
    int ref;
    double elapsed_time;

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    n = id + 1;

    MPI_Reduce(&n, &sum_par, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();

    // printf("n = %d\n", n);
    fflush(stdout);
    MPI_Finalize();
    if (id == 0)
    {
        printf("sum = %d\n", sum_par);

        ref = p * (p + 1) / 2;
        printf("ref = %d\n", ref);
        printf("Elapsed time : %lf\n", elapsed_time);
    }
    return 0;
}