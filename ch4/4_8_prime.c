/*
 * Determine, for all integers less than 1,000,000, 
 * the number of times that two consecutive odd integers are both prime. 
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>

int isPrime(int);

int main(int argc, char *argv[])
{
    int id; /* Process rank */
    int p;  /* Number of processes */
    int cnt = 0;
    int global_cnt = 0;
    int range = 1000000;
    int start; // start number of partition
    int end;   // end number of partition
    double elapsed_time;

    int cur = 0;
    int nex = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    start = (id * range / p) + 1;
    end = (id + 1) * range / p;

    // eliminite even number
    if (start % 2 == 0)
        start++;

    if (id != (p - 1))
        end += 2;

    for (int i = start; i <= end; i += 2)
    {
        // if (isPrime(i) && isPrime(i + 1))
        //     cnt++;
        nex = isPrime(i);
        cnt += cur && nex;
        cur = nex;
        // update next prime
    }

    MPI_Reduce(&cnt, &global_cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    elapsed_time += MPI_Wtime();
    MPI_Finalize();
    if (id == 0)
    {
        printf("%d\n", global_cnt);
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

    for (int i = 5; i <= sqrt_n; i += 6)
    {
        if (n % i == 0 || n % (i + 2) == 0)
            return 0;
    }

    return 1;
}
