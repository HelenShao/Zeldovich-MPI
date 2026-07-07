#pragma once // include this header at most once per file

#include <mpi.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Embedded IC stage API for Abacus (include only from multistep.cpp).
 * @param from_abacus_host Non-zero when MPI is owned by Abacus (skip FFTW teardown in IC driver).
 * @param abacus_comm_2d Abacus 2D cart comm when embedded; MPI_COMM_NULL for standalone.
 */
void IC_InitStage(int from_abacus_host, MPI_Comm abacus_comm_2d);

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
 * @param wisdom_save_dir Abacus ICWisdomSaveDirectory, or NULL for broadcast-only preflight.
 *   Not read from the param header: ICWisdomSaveDirectory is an Abacus parameter (not ZD_*).
 */
int IC_ParamBuffer(const char *bytes, size_t len, const char *wisdom_save_dir);

/** Run full IC generation: argv is Zeldovich_MPI CLI (param_file only; N is derived from NP->ppd). */
int IC_Run(int argc, char **argv);

/** Barriers / hooks after IC_Run (MPI and FFTW state per host contract). */
void IC_FinalizeStage(void);

#ifdef __cplusplus
}
#endif
