/*
 * Write a parallel program to compute the value of pi
 *  using Simpson's Rule
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>

#define n 50

double f(int i)
{
    double x;
    x = (double)i / (double)n;
    return 4.0 / (1.0 + x * x);
}

int main(int argc, char *argv[])
{
    double elapsed_time;
    double area = 0;
    double sum = 0;
    int id;
    int p;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    // sum from n=1 to n/2
    for (int i = id + 1; i <= n / 2; i += p)
        sum += 4.0 * f(2 * i - 1) + 2 * f(2 * i);

    MPI_Reduce(&sum, &area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();
    MPI_Finalize();


    if (id == 0)
    {
        area += f(0) - f(n);
        area /= (3.0 * n);
        printf("Approximation of pi : %13.11lf\n", area);
        printf("Elapsed time : %2.8lf\n", elapsed_time);
    }
    return 0;
}