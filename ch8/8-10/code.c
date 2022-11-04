/**
 * Chapter 8, 8-10
 * Matrix-Vector Multiplication
 * checkerboard biock decomposition of the matrix.
 * read the matrix and the vector from an input file and 
 * print the answer to standard output.
 * 
 * p should be a square number.
 */

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) \
    (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define OPEN_FILE_ERROR -1
#define MALLOC_ERROR -2
#define TYPE_ERROR -3
#define DEBUG 0

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef double dtype;
#define mpitype MPI_DOUBLE

enum MSG_TAG {
    MATRIX_BUFFER,
    BLOCK_VEC,
    PROMPT_MSG
};

int check_square(int p);
void *my_malloc(int id, int bytes);
void read_checkerboard_matrix(char *path, dtype ***mat, int *m, int *n, MPI_Comm grid_comm);
void read_block_vector(char *path, dtype **vec, int *n, MPI_Comm comm);
void print_block_vector(dtype *vec, int n, MPI_Comm comm);
void print_subvector(dtype *vec, int n);
void print_checkerboard_matrix(
    dtype **a,           /* IN -2D matrix */
    int m,               /* IN -Matrix rows */
    int n,               /* IN -Matrix columns */
    MPI_Comm grid_comm); /* IN - Communicator */
void print_submatrix(
    dtype **a, /* OUT - Doubly-subscripted array */
    int rows,  /* OUT - Matrix rows */
    int cols); /* OUT - Matrix cols */

int main(int argc, char *argv[]) {
    double elapsed_time, max_elapsed_time;
    int id, p;
    MPI_Comm grid_comm;           // grid communicator
    MPI_Comm row_comm;            // row communicator
    MPI_Comm col_comm;            // col communicator
    int grid_shape[2] = {0, 0};   // the shape of the grid
    int periodic[2] = {0, 0};     // periodic
    int grid_id;                  // process id in the grid
    int grid_coords[2];           // process's coordinate (grid_coords[0], grid_coords[1])
    dtype **mat;                  // matrix on every process
    int m, n;                     // matrix of m * n
    int local_rows, local_cols;   // process local rows and cols
    dtype *vec;                   // vector on every process
    dtype *local_vec;             // vector on each process
    int vecN;                     // size of the vector
    dtype *local_partial_result;  // local partial result
    dtype *local_result;          // local result
    dtype send_cnt;               // send count
    dtype recv_cnt;               // recv count
    int dest;                     // send destination
    int src;                      // recv source
    int coords[2];
    MPI_Status status;

    if (argc != 4) {
        printf("Usage: ./code <path_to_matrix> <path_to_vector> <path_to_output>\n");
        exit(1);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int is_square = check_square(p);

    // Create the 2D dimension grid, for example:
    // grid_shape = {3, 3}
    MPI_Dims_create(p, 2, grid_shape);
    // Create the 2D Communicator, for example:
    // grid_comm
    // | 0 | 1 | 2 |
    // | 3 | 4 | 5 |
    // | 6 | 7 | 8 |
    MPI_Cart_create(MPI_COMM_WORLD, 2, grid_shape, periodic, 1, &grid_comm);
    // Get the process's id in the grid
    MPI_Comm_rank(grid_comm, &grid_id);
    // Get the process's coordinate in the grid
    MPI_Cart_coords(grid_comm, grid_id, 2, grid_coords);

    if (grid_id == 0) {
        printf("grid_shape = (%d, %d)\n", grid_shape[0], grid_shape[1]);
        fflush(stdout);
    }

    // Get the row communicator (old_comm, color, key, new_comm)
    // color: Process with the same color are in the same new communicator
    // key: process new rank in the communicator
    MPI_Comm_split(grid_comm, grid_coords[0], grid_coords[1], &row_comm);
    // Get the col communicator
    MPI_Comm_split(grid_comm, grid_coords[1], grid_coords[0], &col_comm);
    // Read the matrix from the file in checkerboard decomposition
    read_checkerboard_matrix(argv[1], &mat, &m, &n, grid_comm);
    if (!grid_id) {
        printf("[Matrix] rows = %d, cols = %d\n", m, n);
        fflush(stdout);
    }
    // print_checkerboard_matrix(mat, m, n, grid_comm);

    local_rows = BLOCK_SIZE(grid_coords[0], grid_shape[0], m);
    local_cols = BLOCK_SIZE(grid_coords[1], grid_shape[1], n);

    // The process in the column 0 read the vector
    if (grid_coords[1] == 0) {
        read_block_vector(argv[2], &vec, &vecN, col_comm);
    }
    if (!grid_id) {
        printf("[Vector] size = %d\n", vecN);
        fflush(stdout);
    }
    if (DEBUG) {
        if (grid_coords[1] == 0) {
            print_block_vector(vec, n, col_comm);
        }
    }

    local_partial_result = (dtype *)malloc(local_rows * sizeof(dtype));
    local_result = (dtype *)malloc(local_rows * sizeof(dtype));
    recv_cnt = BLOCK_SIZE(grid_coords[1], grid_shape[1], n);  // id, p, n
    local_vec = (dtype *)malloc(recv_cnt * sizeof(dtype));

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (is_square) {
        // Proc at (i,0) sends subvector to proc at (0,i)
        // Proc at (0,0) just does a copy.
        if (grid_coords[0] == 0 && grid_coords[1] == 0) {  // (0, 0)
            for (int i = 0; i < recv_cnt; ++i) local_vec[i] = vec[i];
        } else if (grid_coords[0] > 0 && grid_coords[1] == 0) {       // (i, 0)
            send_cnt = BLOCK_SIZE(grid_coords[0], grid_shape[0], n);  // send the block vector
            coords[0] = 0;
            coords[1] = grid_coords[0];  // send to (0, i)
            MPI_Cart_rank(grid_comm, coords, &dest);
            MPI_Send(vec, send_cnt, mpitype, dest, BLOCK_VEC, grid_comm);
        } else if (grid_coords[0] == 0 && grid_coords[1] > 0) {  // (0, i)
            coords[0] = grid_coords[1];
            coords[1] = 0;  // recv from (i, 0)
            MPI_Cart_rank(grid_comm, coords, &src);
            MPI_Recv(local_vec, recv_cnt, mpitype, src, BLOCK_VEC, grid_comm, &status);
        }
    } else {
        // Process (0, 0) gathers vector from processes in column 0
        dtype *tmp = (dtype *)malloc(n * sizeof(dtype));
        if (grid_coords[1] == 0) {
            int *gather_recv_cnt, *gather_recv_disp;
            gather_recv_cnt = (int *)malloc(grid_shape[0] * sizeof(int));
            gather_recv_disp = (int *)malloc(grid_shape[0] * sizeof(int));
            for (int i = 0; i < grid_shape[0]; ++i) {
                gather_recv_cnt[i] = BLOCK_SIZE(i, grid_shape[0], n);
            }
            gather_recv_disp[0] = 0;
            for (int i = 1; i < grid_shape[0]; ++i) {
                gather_recv_disp[i] = gather_recv_disp[i - 1] + gather_recv_cnt[i - 1];
            }
            MPI_Gatherv(vec, BLOCK_SIZE(grid_coords[0], grid_shape[0], n), mpitype, tmp, gather_recv_cnt, gather_recv_disp, mpitype, 0, col_comm);
        }

        // Process at (0, 0) scatters vector to row 0
        if (grid_coords[0] == 0) {
            if (grid_shape[1] > 1) {
                int *scatter_send_cnt = (int *)malloc(grid_shape[1] * sizeof(int));
                int *scatter_send_disp = (int *)malloc(grid_shape[1] * sizeof(int));

                for (int i = 0; i < grid_shape[1]; ++i) {
                    scatter_send_cnt[i] = BLOCK_SIZE(i, grid_shape[1], n);
                }
                scatter_send_disp[0] = 0;
                for (int i = 1; i < grid_shape[1]; ++i) {
                    scatter_send_disp[i] = scatter_send_disp[i - 1] + scatter_send_cnt[i - 1];
                }
                MPI_Scatterv(tmp, scatter_send_cnt, scatter_send_disp, mpitype, local_vec, recv_cnt, mpitype, 0, row_comm);

            } else {
                // grid_shape is (i, 0)
                for (int i = 0; i < n; ++i) {
                    local_vec[i] = tmp[i];
                }
            }
        }
    }

    // Rows 0 procs broadcast their local vector to the same column procs
    MPI_Bcast(local_vec, recv_cnt, mpitype, 0, col_comm);

    for (int i = 0; i < local_rows; ++i) {
        local_partial_result[i] = 0.0;
        for (int j = 0; j < local_cols; ++j) {
            local_partial_result[i] += mat[i][j] * local_vec[j];
        }
    }
    MPI_Reduce(local_partial_result, local_result, local_rows, mpitype, MPI_SUM, 0, row_comm);

    // Print block vector
    if (DEBUG) {
        if (!grid_id) {
            printf("[Matrix * Vector] size = %d\n", vecN);
            fflush(stdout);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    MPI_Reduce(&elapsed_time, &max_elapsed_time, 1, mpitype, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!id) {
        printf("(%d * %d), Processes = %d, Time = %12.6f sec\n", m, n, p, max_elapsed_time);
    }

    // write file
    if (grid_coords[1] == 0) {
        // print_block_vector(local_result, n, col_comm);
        write_block_vector(argv[3], local_result, n, col_comm);
    }

    MPI_Finalize();
    return 0;
}

int check_square(int p) {
    int i = 1;
    int is_square = 0;
    while ((i * i) < p) ++i;
    if (i * i == p)
        return 1;
    else
        return 0;
}

/*
 *   Function 'my_malloc' is called when a process wants
 *   to allocate some space from the heap. If the memory
 *   allocation fails, the process prints an error message
 *   and then aborts execution of the program.
 */
void *my_malloc(
    int id,    /* IN - Process rank */
    int bytes) /* IN - Bytes to allocate */
{
    void *buffer;
    if ((buffer = malloc((size_t)bytes)) == NULL) {
        printf("Error: Malloc failed for process %d\n", id);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, MALLOC_ERROR);
    }
    return buffer;
}

void read_checkerboard_matrix(char *path, dtype ***mat, int *m, int *n, MPI_Comm grid_comm) {
    dtype *buffer;  // Read file buffer
    int coords[2];  // Coords of the proc
    int dest;
    int grid_coords[2];
    int grid_id;
    int grid_period[2];
    int grid_shape[2];
    int local_rows;
    int local_cols;
    int p;
    MPI_Status status;
    FILE *infileptr;

    MPI_Comm_rank(grid_comm, &grid_id);
    MPI_Comm_size(grid_comm, &p);

    // Process 0 open files
    // get rows and cols and broadcast
    if (grid_id == 0) {
        infileptr = fopen(path, "r");
        if (infileptr == NULL) {
            printf("Input file '%s' cannot be opened\n", path);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);
        } else {
            fread(m, sizeof(int), 1, infileptr);
            fread(n, sizeof(int), 1, infileptr);
        }
    }
    MPI_Bcast(m, 1, MPI_INT, 0, grid_comm);
    MPI_Bcast(n, 1, MPI_INT, 0, grid_comm);

    MPI_Cart_get(grid_comm, 2, grid_shape, grid_period, grid_coords);

    // Each process determines the size of their local matrix
    local_rows = BLOCK_SIZE(grid_coords[0], grid_shape[0], *m);
    local_cols = BLOCK_SIZE(grid_coords[1], grid_shape[1], *n);

    *mat = (dtype **)my_malloc(grid_id, local_rows * sizeof(dtype *));
    for (int i = 0; i < local_rows; ++i) {
        (*mat)[i] = (dtype *)my_malloc(grid_id, local_cols * sizeof(dtype));
    }

    if (grid_id == 0) {
        buffer = my_malloc(grid_id, *n * sizeof(dtype));
    }

    for (int i = 0; i < grid_shape[0]; ++i) {
        coords[0] = i;
        for (int j = 0; j < BLOCK_SIZE(i, grid_shape[0], *m); ++j) {
            // Read a row in the matrix
            if (grid_id == 0) {
                fread(buffer, sizeof(dtype), *n, infileptr);
            }

            // Distribute the buffer among the process in the row
            for (int k = 0; k < grid_shape[1]; ++k) {
                coords[1] = k;
                dtype *raddr = buffer + BLOCK_LOW(k, grid_shape[1], *n);

                // Get the destination of the process
                MPI_Cart_rank(grid_comm, coords, &dest);

                if (grid_id == 0) {
                    if (dest == 0) {
                        for (int w = 0; w < local_cols; ++w) {
                            (*mat)[j][w] = buffer[w];
                        }
                    } else {
                        // send to dest
                        MPI_Send(raddr, BLOCK_SIZE(k, grid_shape[1], *n), mpitype, dest, MATRIX_BUFFER, grid_comm);
                    }

                } else if (grid_id == dest) {
                    MPI_Recv((*mat)[j], local_cols, mpitype, 0, MATRIX_BUFFER, grid_comm, &status);
                }
            }
        }
    }
    if (grid_id == 0) free(buffer);
}

void read_block_vector(char *path, dtype **vec, int *n, MPI_Comm comm) {
    FILE *infileptr;
    int local_els;
    MPI_Status status;
    int id;
    int p;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &id);

    // Process 0 opens file, broadcasts this value to the other processes
    if (id == p - 1) {
        infileptr = fopen(path, "r");
        if (infileptr == NULL) {
            *n = 0;
        } else {
            fread(n, sizeof(int), 1, infileptr);
        }
    }
    MPI_Bcast(n, 1, MPI_INT, p - 1, comm);
    if (!*n) {
        printf("Input file '%s' cannot be opened\n", path);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);
    }

    local_els = BLOCK_SIZE(id, p, *n);

    // printf("[read_block_vector] vector n = %d\n", *n);
    // fflush(stdout);

    *vec = (dtype *)my_malloc(id, local_els * sizeof(dtype));
    if (id == p - 1) {
        for (int i = 0; i < p - 1; ++i) {
            fread(*vec, sizeof(dtype), BLOCK_SIZE(i, p, *n), infileptr);
            // for (int k = 0; k < BLOCK_SIZE(i, p, *n); ++k) {
            //     printf("[read_block_vector] value n = %6.3f\n", (*vec)[k]);
            //     fflush(stdout);
            // }

            MPI_Send(*vec, BLOCK_SIZE(i, p, *n), mpitype, i, BLOCK_VEC, comm);
        }
        fread(*vec, sizeof(double), BLOCK_SIZE(id, p, *n), infileptr);  // the last read block is it's vec
        fclose(infileptr);
    } else {
        MPI_Recv(*vec, BLOCK_SIZE(id, p, *n), mpitype, p - 1, BLOCK_VEC, comm, &status);
    }
}

void print_block_vector(dtype *vec, int n, MPI_Comm comm) {
    int id;      // process rank
    int p;       // processes
    int prompt;  // Dummy
    dtype *tmp;  // Other process's subvector
    MPI_Status status;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &id);

    if (!id) {
        print_subvector(vec, BLOCK_SIZE(id, p, n));
        if (p > 1) {
            tmp = (dtype *)my_malloc(id, BLOCK_SIZE(p - 1, p, n) * sizeof(dtype));
            for (int i = 1; i < p; ++i) {
                MPI_Send(&prompt, 1, MPI_INT, i, PROMPT_MSG, comm);
                MPI_Recv(tmp, BLOCK_SIZE(i, p, n), mpitype, i, BLOCK_VEC, comm, &status);
                print_subvector(tmp, BLOCK_SIZE(i, p, n));
            }
            free(tmp);
        }
        printf("\n");
    } else {
        MPI_Recv(&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm, &status);
        MPI_Send(vec, BLOCK_SIZE(id, p, n), mpitype, 0, BLOCK_VEC, comm);
    }
}

void print_subvector(dtype *vec, int n) {
    for (int i = 0; i < n; ++i) {
        printf("%.8f ", vec[i]);
    }
}

void write_block_vector(char *path, dtype *vec, int n, MPI_Comm comm) {
    int id;      // process rank
    int p;       // processes
    int prompt;  // Dummy
    dtype *tmp;  // Other process's subvector
    MPI_Status status;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &id);

    if (!id) {
        FILE *foutptr = fopen(path, "w");
        fprintf(foutptr, "%d\n", n);
        write_subvector(foutptr, vec, BLOCK_SIZE(id, p, n));
        if (p > 1) {
            tmp = (dtype *)my_malloc(id, BLOCK_SIZE(p - 1, p, n) * sizeof(dtype));
            for (int i = 1; i < p; ++i) {
                MPI_Send(&prompt, 1, MPI_INT, i, PROMPT_MSG, comm);
                MPI_Recv(tmp, BLOCK_SIZE(i, p, n), mpitype, i, BLOCK_VEC, comm, &status);
                write_subvector(foutptr, tmp, BLOCK_SIZE(i, p, n));
            }
            free(tmp);
        }
        fflush(foutptr);
        fclose(foutptr);
    } else {
        MPI_Recv(&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm, &status);
        MPI_Send(vec, BLOCK_SIZE(id, p, n), mpitype, 0, BLOCK_VEC, comm);
    }
}
void write_subvector(FILE *fp, dtype *vec, int n) {
    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%.8f ", vec[i]);
    }
}

void print_checkerboard_matrix(dtype **a, int m, int n, MPI_Comm grid_comm) {
    void *buffer;       /* Room to hold 1 matrix row */
    int coords[2];      /* Grid coords of process
                                 sending elements */
    int datum_size;     /* Bytes per matrix element */
    int els;            /* Elements received */
    int grid_coords[2]; /* Coords of this process */
    int grid_id;        /* Process rank in grid */
    int grid_period[2]; /* Wraparound */
    int grid_size[2];   /* Dims of process grid */
    int i, j, k;
    void *laddr;       /* Where to put subrow */
    int local_cols;    /* Matrix cols on this proc */
    int p;             /* Number of processes */
    int src;           /* ID of proc with subrow */
    MPI_Status status; /* Result of receive */

    MPI_Comm_rank(grid_comm, &grid_id);
    MPI_Comm_size(grid_comm, &p);
    datum_size = sizeof(dtype);

    MPI_Cart_get(grid_comm, 2, grid_size, grid_period,
                 grid_coords);
    local_cols = BLOCK_SIZE(grid_coords[1], grid_size[1], n);

    if (!grid_id)
        buffer = my_malloc(grid_id, n * datum_size);

    /* For each row of the process grid */
    for (i = 0; i < grid_size[0]; i++) {
        coords[0] = i;

        /* For each matrix row controlled by the process row */
        for (j = 0; j < BLOCK_SIZE(i, grid_size[0], m); j++) {
            /* Collect the matrix row on grid process 0 and
            print it */
            if (!grid_id) {
                for (k = 0; k < grid_size[1]; k++) {
                    coords[1] = k;
                    MPI_Cart_rank(grid_comm, coords, &src);
                    els = BLOCK_SIZE(k, grid_size[1], n);
                    laddr = buffer +
                            BLOCK_LOW(k, grid_size[1], n) * datum_size;
                    if (src == 0) {
                        memcpy(laddr, a[j], els * datum_size);
                    } else {
                        MPI_Recv(laddr, els, mpitype, src, 0,
                                 grid_comm, &status);
                    }
                }
                print_subvector(buffer, n);
                putchar('\n');
            } else if (grid_coords[0] == i) {
                MPI_Send(a[j], local_cols, mpitype, 0, 0,
                         grid_comm);
            }
        }
    }
    if (!grid_id) {
        free(buffer);
    }
}

void print_submatrix(dtype **a, int rows, int cols) {
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%6.3f ", a[i][j]);
        }
        putchar('\n');
    }
}
