#ifndef FFT_WISDOM_H
#define FFT_WISDOM_H

#include <mpi.h>
#include "../precision.h"

// Default wisdom file (relative to cwd) if FFTW_WISDOM_FILE is unset.
#define FFTW_WISDOM_FILENAME "fftw_wisdom_float"

#ifdef __cplusplus
extern "C" {
#endif

// Import: rank 0 reads FFTW_WISDOM_FILE if set (absolute path from job script),
// else FFTW_WISDOM_FILENAME; then broadcast wisdom to all ranks in comm.
void fft_wisdom_import_broadcast(int rank, MPI_Comm comm);

// Export: rank 0 writes to the same path. Rank 0 prints a short wisdom preview first.
void fft_wisdom_export_rank0(int rank);

#ifdef __cplusplus
}
#endif

#endif

