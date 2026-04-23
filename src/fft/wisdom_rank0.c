/*
 * Rank-0 FFTW planning and wisdom export
 */

#include "wisdom_rank0.h"
#include "fft_wisdom.h"
#include "../config.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// IMPORT WISDOM FROM FILE
static void wisdom_import_file(void)
{
    const int imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(FFTW_WISDOM_FILENAME);
    if (!imported) {
        FFTW_FORGET_WISDOM();
        printf("[wisdom_rank0] No usable wisdom at '%s' (%s); measuring cold.\n",
               FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stdout);
    } else {
        printf("[wisdom_rank0] Imported wisdom from '%s' (%s precision).\n", FFTW_WISDOM_FILENAME,
               PRECISION_NAME);
        fflush(stdout);
    }
}

// EXPORT WISDOM TO FILE
static int wisdom_export_file(void)
{
    const int ok = FFTW_EXPORT_WISDOM_TO_FILENAME(FFTW_WISDOM_FILENAME);
    if (!ok) {
        fprintf(stderr, "[wisdom_rank0] Failed to export wisdom to '%s' (%s precision)\n",
                FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stderr);
        return -1;
    }
    printf("[wisdom_rank0] Exported wisdom to '%s' (%s precision)\n", FFTW_WISDOM_FILENAME,
           PRECISION_NAME);
    fflush(stdout);
    return 0;
}

// PLAN AND EXPORT WISDOM
int wisdom_rank0_plans_and_export(int N, int narray, fftw_complex_t *plan_buffer,
                                  fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *dummy_1d = NULL; // dummy 1D buffer for 1D FFT plan

    if (plan_2d_out == NULL || plan_1d_out == NULL) {
        fprintf(stderr, "[wisdom_rank0] output plan pointers must be non-NULL\n");
        return -1;
    }
    *plan_2d_out = NULL;
    *plan_1d_out = NULL;

    if (N <= 0 || (N % 4) != 0) {
        fprintf(stderr,
                "[wisdom_rank0] N (PPD) must be divisible by 4 (got N=%d)\n",
                N);
        return -1;
    }
    if (plan_buffer == NULL) {
        fprintf(stderr, "[wisdom_rank0] plan_buffer must be non-NULL\n");
        return -1;
    }

    /* Keep same order as fft_setup.c: thread init before wisdom import. */
    static int fftw_threads_initialized = 0;
    int fft_threads_2d = omp_get_max_threads();
    if (!fftw_threads_initialized) {
        fft_threads_2d = omp_get_max_threads();

        if (FFTW_INIT_THREADS() == 0) {
            fprintf(stderr, "[wisdom_rank0] FFTW_INIT_THREADS failed (%s precision)\n",
                    PRECISION_NAME);
            return -1;
        }
        FFTW_PLAN_WITH_NTHREADS(fft_threads_2d);
        printf("[wisdom_rank0] %s precision: 2D plan uses %d FFTW threads; 1D plan will use 1 thread\n",
               PRECISION_NAME, fft_threads_2d);
        fflush(stdout);
        fftw_threads_initialized = 1;
    }

    wisdom_import_file();

    {
        int n[2] = { N, N };
        *plan_2d_out = FFTW_PLAN_MANY_DFT(
            2,
            n,
            narray,
            plan_buffer,
            NULL,
            1,
            N * N,
            plan_buffer,
            NULL,
            1,
            N * N,
            FFT_SIGN,
            FFTW_MEASURE);
    }

    printf("[wisdom_rank0] 2D planner thread target (OMP max): %d\n", omp_get_max_threads());
    fflush(stdout);

    // Force 1 thread for 1D plan creation so wisdom matches runtime policy.
    FFTW_PLAN_WITH_NTHREADS(1);

    if (posix_memalign((void **)&dummy_1d, ALIGN_BYTES, sizeof(fftw_complex_t) * (size_t)N) != 0) {
        fprintf(stderr, "[wisdom_rank0] posix_memalign failed for 1D dummy\n");
        if (*plan_2d_out) {
            FFTW_DESTROY_PLAN(*plan_2d_out);
            *plan_2d_out = NULL;
        }
        return -1;
    }

    *plan_1d_out =
        FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, FFT_SIGN, FFTW_MEASURE);

    printf("[wisdom_rank0] 1D plan created with FFTW threads: 1\n");
    fflush(stdout);
    free(dummy_1d);

    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[wisdom_rank0] plan creation failed (NULL plan)\n");
        if (*plan_2d_out) {
            FFTW_DESTROY_PLAN(*plan_2d_out);
            *plan_2d_out = NULL;
        }
        if (*plan_1d_out) {
            FFTW_DESTROY_PLAN(*plan_1d_out);
            *plan_1d_out = NULL;
        }
        return -1;
    }

    if (wisdom_export_file() != 0) {
        FFTW_DESTROY_PLAN(*plan_2d_out);
        FFTW_DESTROY_PLAN(*plan_1d_out);
        *plan_2d_out = NULL;
        *plan_1d_out = NULL;
        return -1;
    }

    return 0;
}
