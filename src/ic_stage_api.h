#pragma once // include this header at most once per file

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Embedded IC stage API for Abacus (include only from multistep.cpp).
 * @param from_abacus_host Non-zero when MPI is owned by Abacus (skip FFTW teardown in IC driver).
 */
void IC_InitStage(int from_abacus_host);

/**
 * Rank-0 FFTW wisdom preflight (standalone CLI): reads param_file on rank 0.
 * Embedded Abacus hosts should use IC_ParamBuffer instead (no second file read).
 */
int IC_Rank0Wisdom(const char *param_file);

/**
 * Embedded IC entry (Abacus): rank-0 param header bytes already broadcast by host.
 *
 * Wisdom preflight (before IC driver), rank 0 only:
 * - Parse params from memory, plan 2D+1D FFTW; optionally write fftw_wisdom.wisdom.
 * - All ranks barrier; wisdom MPI_Bcast happens in setup_fftw_plans_full (driver).
 *
 * @param param_path Path string for relative par2 lookups only (not re-read on rank 0).
 * @param wisdom_save_dir Abacus ICWisdomSaveDirectory, or NULL for broadcast-only preflight.
 */
int IC_ParamBuffer(const char *bytes, size_t len, const char *param_path, const char *wisdom_save_dir);

/** Run full IC generation: argv is Zeldovich_MPI CLI (param_file only; N is derived from NP->ppd). */
int IC_Run(int argc, char **argv);

/** Barriers / hooks after IC_Run (MPI and FFTW state per host contract). */
void IC_FinalizeStage(void);

#ifdef __cplusplus
}
#endif
