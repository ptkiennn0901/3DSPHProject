/* hello_mpi.c — Test MPI co ban. Moi rank in rank cua no va hostname cua may.
 *
 *   mpicc -O2 -o hello src/hello_mpi.c
 *   mpirun --hostfile ~/hostfile -np 6 ./hello */
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    int rank;
    int size;
    char host[256];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    gethostname(host, sizeof(host));
    printf("Xin chao tu rank %d / %d tren may '%s'\n", rank, size, host);

    MPI_Finalize();
    return 0;
}
