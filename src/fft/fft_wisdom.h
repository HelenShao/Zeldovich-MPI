#ifndef FFT_WISDOM_H
#define FFT_WISDOM_H

#include "../precision.h"
#include <mpi.h>

// Wisdom file path (relative to cwd), same for wisdom_rank0 / IC_Rank0Wisdom and IC_Run import.
#define FFTW_WISDOM_FILENAME "fftw_wisdom.wisdom"

#ifdef __cplusplus
extern "C" {
#endif

// Rank 0 imports FFTW_WISDOM_FILENAME (from preflight), broadcasts the wisdom string,
// and each rank loads it via FFTW_IMPORT_WISDOM_FROM_STRING.
// local_wisdom_dir is unused (kept for API compatibility).
// Returns 0 on success, non-zero on failure.
int fft_wisdom_import_rank0_broadcast_local(int rank, MPI_Comm comm, const char *local_wisdom_dir);

#ifdef __cplusplus
}
#endif

#endif
