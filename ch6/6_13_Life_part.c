/*
 * Life is an example of a cellular automaton. It consists of a rectangular grid of cells
 * During each iteration,
 *
 * 1. A dead cell with exactly three neighbors becomes a live cell.
 * 2. A live cell with two or three neighbors stays alive.
 * 3. A live cell with less than two neighbors or more than three neighbors becomes a dead cell.
 *
 * Write a parallel program that reads from a file an m  x  n matrix containing
 * the initial state of the game. It should play the game of Life for j iterations,
 * printing the state of the game once every k iterations, where j and k are command-line
 * arguments.
 *
 *  ./Life <j> <k>
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>

#define LIVE 1
#define DEAD 0

#define MAX_N 10

#define M 5
#define N 5
#define TEST_MATRIX { \
    {1, 0, 0, 1, 1},  \
    {0, 0, 0, 1, 1},  \
    {0, 0, 0, 0, 1},  \
    {1, 1, 1, 0, 1},  \
    {0, 0, 0, 0, 1}};

void print_state(int, int, int **);

int count_live(int **, int, int);

int main(int argc, char *argv[])
{
    int id; /* Process rank */
    int p;  /* Number of processes */
    int j;
    int k;

    int state_matrix[MAX_N][MAX_N] = TEST_MATRIX;
    int temp_matrix[MAX_N][MAX_N];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    j = atoi(argv[1]);
    k = atoi(argv[2]);

    // consider i-1, i, i+1 row
    // start iteration
    for (int i = id; i < j; j += p)
    {
        // case 1
        if (count_live(state_matrix, i, j) < 2)
        {
            temp_matrix[i][j]
        }
    }

    // boardcast every iteration

    MPI_Reduce(&cnt, &global_cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Finalize();

    if (id == 0)
    {
        printf("%d\n", global_cnt);
    }
    return 0;
}

void print_state(int m, int n, int **matrix)
{
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (matrix[i][j] != 0)
                printf("O ");
            else
                printf("X ");
        }
        printf("\n");
    }
}

int count_live(int **matrix, int i, int j)
{
    int cnt = 0;
    cnt += (matrix[i - 1][j] == LIVE) ? 1 : 0;
    cnt += (matrix[i][j - 1] == LIVE) ? 1 : 0;
    cnt += (matrix[i][j + 1] == LIVE) ? 1 : 0;
    cnt += (matrix[i + 1][j] == LIVE) ? 1 : 0;
    return cnt;
}