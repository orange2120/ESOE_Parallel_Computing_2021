/*
 * 
 *  ./Life <filename> <j> <k>
 */

#define DATA_MSG 0

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) \
    (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

int main(int argc, char *argv[])
{
    int **a;
    int *storage; // local portion of array elements
    int i, j, k;
    int iterations, printFreq;
    int id;
    int m, n;
    int p;
    double elapsed_time;
    MPI_Status status; /* Result of receive */

    FILE *fp; /* Input file pointer */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (id == p - 1)
    {
        // read m, n;
        fp = fopen(argv[1], "r");
        fscanf(fp, "%d", &m);
        fscanf(fp, "%d", &n);

        // read iteration, printFreq
        iterations = atoi(argv[2]);
        printFreq = atoi(argv[3]);
    }

    MPI_Bcast(&m, 1, MPI_INT, p - 1, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, p - 1, MPI_COMM_WORLD);
    MPI_Bcast(&iterations, 1, MPI_INT, p - 1, MPI_COMM_WORLD);
    MPI_Bcast(&printFreq, 1, MPI_INT, p - 1, MPI_COMM_WORLD);

    int **matrix; // global
    int localRowStart = BLOCK_LOW(id, p, m);
    int localRowEnd = BLOCK_HIGH(id, p, m);
    int localRowCnt = BLOCK_SIZE(id, p, m);
    // allocate more two lines for overlap (calculate 8 neigbors)
    int **localRows = (int **)malloc((localRowCnt + 2) * sizeof(int *));
    for (i = 0; i < localRowCnt + 2; ++i)
        localRows[i] = (int *)malloc(n * sizeof(int));

    if (id == p - 1)
    {
        // read whole matrix for last process
        matrix = (int **)malloc(m * sizeof(int *));
        for (i = 0; i < m; ++i)
        {
            matrix[i] = (int *)malloc(n * sizeof(int *));
            for (j = 0; j < n; ++j)
            {
                fscanf(fp, "%d", &matrix[i][j]);
            }
        }
        fclose(fp);
    }

    for (k = 1; k <= iterations; ++k)
    {
        // send to process i
        if (id == p - 1)
        {
            // send more two lines to process i
            for (i = 0; i < p - 1; ++i)
            {
                int start = BLOCK_LOW(i, p, m);
                int end = BLOCK_HIGH(i, p, m);
                for (j = start - 1; j <= end + 1; ++j)
                {
                    if (j < 0 || j >= m)
                        continue;
                    MPI_Send(matrix[j], n, MPI_INT, i, DATA_MSG, MPI_COMM_WORLD);
                }
            }

            // for last process
            for (i = localRowStart - 1; i <= localRowEnd; ++i)
            {
                for (j = 0; j < n; ++j)
                    localRows[i - localRowStart + 1][j] = matrix[i][j];
            }
        }
        else
        {
            // recv from process p - 1
            for (i = localRowStart - 1; i <= localRowEnd + 1; ++i)
            {
                if (i < 0 || i >= m)
                    continue;
                MPI_Recv(localRows[i - localRowStart + 1], n, MPI_INT, p - 1, DATA_MSG, MPI_COMM_WORLD, &status);
            }
        }

        // set margin in process 0 and process p - 1 to -1
        if (id == 0)
        {
            for (j = 0; j < n; ++j)
                localRows[0][j] = -1; // first line
        }
        if (id == p - 1)
        {
            for (j = 0; j < n; ++j)
                localRows[localRowCnt + 1][j] = -1; // last line
        }

        // printf("Process %d recv %d rows: \n", id, localRowCnt + 2);
        // for (i = 0; i < localRowCnt + 2; ++i) {
        //   for (j = 0; j < n; ++j) {
        //     printf("%d ", localRows[i][j]);
        //   }
        //   printf("\n");
        // }
        // fflush(stdout);

        // do calculation
        int **nextState = (int **)malloc(localRowCnt * sizeof(int *));
        for (i = 1; i < localRowCnt + 1; ++i)
        {
            nextState[i - 1] = (int *)malloc(n * sizeof(int));
            for (j = 0; j < n; ++j)
            {
                int x, y;
                int life = 0;
                for (x = i - 1; x <= i + 1; ++x)
                {
                    for (y = j - 1; y <= j + 1; ++y)
                    {
                        if (y < 0 || y >= n)
                            continue;
                        if (x == i && y == j)
                            continue;
                        if (localRows[x][y] == 1)
                            ++life;
                    }
                }
                if (localRows[i][j] == 0 && life == 3)
                    nextState[i - 1][j] = 1;
                else if (localRows[i][j] == 1 && (life == 2 || life == 3))
                    nextState[i - 1][j] = 1;
                else
                    nextState[i - 1][j] = 0;
            }
        }

        // printf("Process %d nextState: \n", id);
        // for (i = 0; i < localRowCnt; ++i) {
        //   for (j = 0; j < n; ++j) {
        //     printf("%d", nextState[i][j]);
        //   }
        //   printf("\n");
        // }
        // printf("\n");
        // fflush(stdout);

        // recv from process i
        if (id == p - 1)
        {
            for (i = 0; i < p - 1; ++i)
            {
                int start = BLOCK_LOW(i, p, m);
                int end = BLOCK_HIGH(i, p, m);
                for (j = start; j <= end; ++j)
                {
                    MPI_Recv(matrix[j], n, MPI_INT, i, DATA_MSG, MPI_COMM_WORLD, &status);
                }
            }

            // for last process, copy nextState to matrix
            for (i = localRowStart; i <= localRowEnd; ++i)
            {
                for (j = 0; j < n; ++j)
                    matrix[i][j] = nextState[i - localRowStart][j];
            }
        }
        else
        {
            // send to process p - 1
            for (i = 0; i < localRowCnt; ++i)
            {
                MPI_Send(nextState[i], n, MPI_INT, p - 1, DATA_MSG, MPI_COMM_WORLD);
            }
        }

        // print answer
        if (id == p - 1 && k % printFreq == 0)
        {
            fprintf(stderr, "After %d iteration:\n", k);
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < n; ++j)
                {
                    fprintf(stderr, "%d", matrix[i][j]);
                }
                fprintf(stderr, "\n");
            }
            fflush(stdout);
        }
    }

    elapsed_time += MPI_Wtime();
    if (id == p - 1)
    {
        printf("Running with process number: %d\n", p);
        printf("Execution time %8.6f\n", elapsed_time);
    }

    MPI_Finalize();
}