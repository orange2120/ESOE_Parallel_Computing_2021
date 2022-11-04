#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)

const int s = 2;
const double d = 0.3;

int main(int argc, char* argv[]) {
    const long N_SAMPLE = atol(argv[1]);
    MPI_Init(&argc, &argv);

    int id, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Barrier(MPI_COMM_WORLD);
    double local_elapsed_time = -MPI_Wtime(), max_elapsed_time;
    long start = BLOCK_LOW(id, p, N_SAMPLE);
    long end = BLOCK_HIGH(id, p, N_SAMPLE);

    long local = 0, total = 0;  // local valid sample count, total valid sample count
    srand((unsigned)time(NULL) + id);

    const double threshold = d / 2;
    for (long i = start; i <= end; ++i) {
        double x = (double)rand() / RAND_MAX * s;
        double y = (double)rand() / RAND_MAX * s;
        double z = (double)rand() / RAND_MAX * s;
        double diag = sqrt(x * x + y * y + z * z);
        if (sin(acos((x + y + z) / (sqrt(3) * diag))) * diag > d) {
            ++local;
        }
    }
    MPI_Reduce(&local, &total, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    double volume = 0;
    if (id == 0) {
        volume = s * s * s * ((double)total / (double)N_SAMPLE);
    }
    local_elapsed_time += MPI_Wtime();
    MPI_Reduce(&local_elapsed_time, &max_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (id == 0) {
        printf(
            "Volume: %.6lf, time elapsed = "
            "%.8lf\n",
            volume, max_elapsed_time);
        fflush(stdout);
    }
    MPI_Finalize();
}
