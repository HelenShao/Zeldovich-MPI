#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** Full Zeldovich_MPI pipeline; requires MPI initialized on MPI_COMM_WORLD. */
int zeldovich_mpi_driver_run(int argc, char **argv);

#ifdef __cplusplus
}
#endif
