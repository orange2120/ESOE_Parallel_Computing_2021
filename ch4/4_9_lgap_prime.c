
/*
 * Determine, for all integers less than 1,000,000,
 * the largest gap between a pair of consecutive prime numbers.
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>

int isPrime(int);

int main(int argc, char *argv[])
{
    int id; /* Process rank */
    int p;  /* Number of processes */
    int gap = 0;
    int global_gap = 0;
    int range = 1000000;
    int start; // start number of partition
    int end;   // end number of partition
    int prev_prime = 0;

    double elapsed_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    start = (id * range / p) + 1;
    end = (id + 1) * range / p;

    // keep start, end are odd
    // eliminite even number
    if (start % 2 == 0)
        start++;

    if (id != (p - 1))
        end += 2;

    // printf("start=%d end = %d\n", start, end);
    for (int i = start; prev_prime <= end && i < range; i += 2)
    {
        if (isPrime(i))
        {
            if (prev_prime)
            {
                if ((i - prev_prime) > gap)
                {
                    gap = i - prev_prime;
                }
            }
            prev_prime = i;
        }
    }

    MPI_Reduce(&gap, &global_gap, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    MPI_Finalize();
    if (id == 0)
    {
        printf("%d\n", global_gap);
        printf("Elapsed time : %lf\n", elapsed_time);
    }
    return 0;
}

int isPrime(int n)
{
    int sqrt_n = sqrt(n);
    if (n <= 1)
        return 0;
    if (n <= 3)
        return 1;
    if (n % 2 == 0 || n % 3 == 0)
        return 0;

    for (int i = 5; i < sqrt_n; i += 6)
    {
        if (n % i == 0 || n % (i + 2) == 0)
            return 0;
    }

    return 1;
}
