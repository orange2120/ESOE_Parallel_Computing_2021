#include <stdio.h>

#include "mpi.h"

int flip(int **matrix, int m, int n, int x, int y)
{
    int life = 0;
    for (int i = x - 1; i <= x + 1; ++i)
    {
        for (int j = y - 1; j <= y + 1; ++j)
        {
            if (i < 0 || i >= m)
                continue;
            if (j < 0 || j >= n)
                continue;
            if (i == x && j == y)
                continue;
            if (matrix[i][j])
                ++life;
        }
    }

    if (!matrix[x][y] && life == 3)
        return 1;
    if (matrix[x][y] && (life == 2 || life == 3))
        return 1;
    return 0;
}

int main(int argc, char *argv[])
{
    int m, n;
    int p, id;
    FILE *fp;
    fp = fopen(argv[1], "r");
    fscanf(fp, "%d", &m);
    fscanf(fp, "%d", &n);
    int iterations = atoi(argv[2]);
    int printFreq = atoi(argv[3]);

    int **matrix = (int **)malloc(m * sizeof(int *));
    for (int i = 0; i < m; ++i)
    {
        matrix[i] = (int *)malloc(n * sizeof(int));
        for (int j = 0; j < n; ++j)
        {
            fscanf(fp, "%d", &matrix[i][j]);
        }
    }
    fclose(fp);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_time = -MPI_Wtime();

    for (int k = 1; k <= iterations; ++k)
    {
        int tmpMatrix[m][n];
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                tmpMatrix[i][j] = flip(matrix, m, n, i, j);
            }
        }

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                matrix[i][j] = tmpMatrix[i][j];
            }
        }
        if (k % printFreq == 0)
        {
            fprintf(stderr, "After %d iteration:\n", k);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    fprintf(stderr, "%d", tmpMatrix[i][j]);
                }
                fprintf(stderr, "\n");
            }
        }
    }
    elapsed_time += MPI_Wtime();
    printf("Running with process nuber: %d, ", p);
    printf("Execution time %8.6f\n", elapsed_time);

    MPI_Finalize();
}