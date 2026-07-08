#include <mpi.h>
#include <stdbool.h>

bool zd_log_on_this_rank(void)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank == 0;
}
