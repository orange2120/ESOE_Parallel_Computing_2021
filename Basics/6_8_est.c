/*
 * Estimate the network latency and bandwidth.
 */

#include "mpi.h"
#include <stdio.h>

#define MESSAGE_SZ 1 // bytes

int main(int argc, char *argv[])
{
    int id; /* Process rank */
    int p;  /* Number of processes */
    double elapsed_time;
    double latency;
    double send_time;
    char msg_tx[MESSAGE_SZ];
    char msg_rx[MESSAGE_SZ];
    MPI_Status status;
    static int cnt = 0;

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (id == 0)
    {
        // send message to proc 1
        MPI_Send(msg_tx, MESSAGE_SZ, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
        send_time = MPI_Wtime();
        MPI_Recv(msg_rx, MESSAGE_SZ, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);

        latency = (MPI_Wtime() - send_time) / 2;
        elapsed_time += MPI_Wtime();
        // cnt++;
    }
    else
    {
        MPI_Recv(msg_rx, MESSAGE_SZ, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Send(msg_rx, MESSAGE_SZ, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    if (id == 0)
    {
        printf("Latency = %lf\n", latency);
        printf("time elapse = %lf\n", elapsed_time);
    }     

    return 0;
}