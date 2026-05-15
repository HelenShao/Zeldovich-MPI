#ifndef FFT_WISDOM_H
#define FFT_WISDOM_H

#include "../precision.h"
#include <mpi.h>

// Wisdom file path (relative to cwd), same for wisdom_rank0 and MPI run.
#define FFTW_WISDOM_FILENAME "fftw_wisdom_float.wisdom"

#ifdef __cplusplus
extern "C" {
#endif

// Rank 0 imports FFTW_WISDOM_FILENAME, broadcasts wisdom string, and each rank
// writes/imports its own local wisdom file in local_wisdom_dir.
// Returns 0 on success, non-zero on failure.
int fft_wisdom_import_rank0_broadcast_local(int rank, MPI_Comm comm, const char *local_wisdom_dir);

#ifdef __cplusplus
}
#endif

#endif
