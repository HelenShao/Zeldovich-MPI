#ifndef FFT_WISDOM_H
#define FFT_WISDOM_H

#include <mpi.h>
#include "../precision.h"

// Shared filename for FFTW wisdom (precision-dependent via FFTW_* macros).
// For single-precision (default), this will store fftwf_* wisdom.
#define FFTW_WISDOM_FILENAME "fftw_wisdom_float"

#ifdef __cplusplus
extern "C" {
#endif

// Import per-rank FFTW wisdom from "<FFTW_WISDOM_FILENAME>_rankXXXX".
// Broadcasting is intentionally disabled; each rank manages its own file.
void fft_wisdom_import_per_rank(int rank, MPI_Comm comm);

// Export per-rank FFTW wisdom to "<FFTW_WISDOM_FILENAME>_rankXXXX".
// This creates the file if it did not exist, and overwrites it otherwise.
void fft_wisdom_export_per_rank(int rank);

#ifdef __cplusplus
}
#endif

#endif

