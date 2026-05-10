#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Embedded IC stage API for Abacus (include only from multistep.cpp).
 * @param from_abacus_host Non-zero when MPI is owned by Abacus (skip FFTW teardown in driver).
 */
void IC_InitStage(int from_abacus_host);

/** Optional rank-0 FFTW wisdom hook (reserved; may be no-op). */
void IC_Rank0Wisdom(void);

/** Run full IC generation: argv is the usual Zeldovich_MPI CLI (N, param_file). */
int IC_Run(int argc, char **argv);

/** Barriers / hooks after IC_Run (MPI and FFTW state per host contract). */
void IC_FinalizeStage(void);

#ifdef __cplusplus
}
#endif
