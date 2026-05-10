#include "ic_stage_api.h"

#include <mpi.h>

#include "ic_embed_flags.h"
#include "zeldovich_mpi_driver.h"

bool zeldovich_ic_embedded = false;

void IC_InitStage(int from_abacus_host)
{
    zeldovich_ic_embedded = (from_abacus_host != 0);
}

void IC_Rank0Wisdom(void)
{
    /* Host may run wisdom_rank0 preflight; reserved for future in-process hook. */
}

int IC_Run(int argc, char **argv)
{
    return zeldovich_mpi_driver_run(argc, argv);
}

void IC_FinalizeStage(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    zeldovich_ic_embedded = false;
}
