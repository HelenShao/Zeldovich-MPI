#include <stdio.h>
#include <mpi.h>

#include "ic_embed_flags.h"
#include "zeldovich_mpi_driver.h"

int main(int argc, char **argv)
{
    int provided = MPI_THREAD_SINGLE;
    const int required = MPI_THREAD_SINGLE;
    const int ret = MPI_Init_thread(nullptr, nullptr, required, &provided);
    if (ret != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init_thread failed with error code %d\n", ret);
        return 1;
    }

    zeldovich_ic_embedded = false;
    const int r = zeldovich_mpi_driver_run(argc, argv);

    MPI_Finalize();
    return r;
}
