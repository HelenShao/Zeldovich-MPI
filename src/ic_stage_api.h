#pragma once // include this header at most once per file

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Embedded IC stage API for Abacus (include only from multistep.cpp).
 * @param from_abacus_host Non-zero when MPI is owned by Abacus (skip FFTW teardown in IC driver).
 */
void IC_InitStage(int from_abacus_host);

/**
 * Rank-0 FFTW wisdom preflight (embedded IC): same work as standalone `wisdom_rank0`.
 * Loads param_file w/ ZD param parser to get PPD and narray (qdensity/qPLT),
 * measures/plans on rank 0 and writes `FFTW_WISDOM_FILENAME` in cwd
 * (must match where fft_wisdom_import_rank0_broadcast_local() will look during IC_Run).
 * Other ranks only participate in MPI_Barrier.
 */
int IC_Rank0Wisdom(const char *param_file);

/** Run full IC generation: argv is Zeldovich_MPI CLI (param_file only; N is derived from NP->ppd). */
int IC_Run(int argc, char **argv);

/** Barriers / hooks after IC_Run (MPI and FFTW state per host contract). */
void IC_FinalizeStage(void);

#ifdef __cplusplus
}
#endif
