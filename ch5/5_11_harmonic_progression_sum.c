/*
 * Write a parallel program that computes sums to arbitrary 
 * precision after the decimal point
 * ./hps <n> <d>
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
    int id; /* Process rank */
    int p;  /* Number of processes */
    int n = 1;
    int d = 1; // digits of decimal

    double sum = 0;
    double global_sum = 0;
    double elapsed_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    n = atoi(argv[1]);
    d = atoi(argv[2]);

    // boardcast n, d to process
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // printf("n = %d, d = %d\n", n, d);

    for (int i = id + 1; i <= n; i += p)
    {
        sum += 1.0 / (double)i;
    }

    MPI_Reduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    MPI_Finalize();

    if (id == 0)
    {
        printf("sum = %.*f\n", d, global_sum);
        printf("Elapsed time : %lf\n", elapsed_time);
    }
    return 0;
}
