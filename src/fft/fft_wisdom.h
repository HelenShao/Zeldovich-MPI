#ifndef FFT_WISDOM_H
#define FFT_WISDOM_H

#include "../precision.h"
#include <mpi.h>
#include <stdbool.h>

// Basename for rank-0 preflight export (directory set via zd_wisdom_set_dir).
#define FFTW_WISDOM_BASENAME "fftw_wisdom.wisdom"
#define FFTW_WISDOM_FILENAME FFTW_WISDOM_BASENAME

#ifdef __cplusplus
extern "C" {
#endif

// Set directory for rank-0 preflight export/import (NULL or "" -> cwd basename only).
void zd_wisdom_set_dir(const char *dir);

// Full path: "{dir}/fftw_wisdom.wisdom" or "fftw_wisdom.wisdom" when dir unset.
const char *zd_wisdom_rank0_file(void);

// Rank 0: mkdir parent of zd_wisdom_rank0_file() when it contains a '/'. Returns 0 on success.
int zd_wisdom_ensure_dir_rank0(void);

// Set after successful MPI wisdom broadcast (driver); skips duplicate broadcast.
void zd_wisdom_set_preflight_broadcast_done(bool done);
bool zd_wisdom_preflight_broadcast_done(void);

// MPI_Bcast wisdom from rank 0; all ranks IMPORT_FROM_STRING. No-op if already broadcast.
// Call from setup_fftw_plans_full after FFTW_INIT_THREADS on all ranks.
// local_wisdom_dir: file fallback on rank 0 when FFTW has no in-memory wisdom (standalone).
int fft_wisdom_broadcast_from_rank0(int rank, MPI_Comm comm, const char *local_wisdom_dir);

#ifdef __cplusplus
}
#endif

#endif
