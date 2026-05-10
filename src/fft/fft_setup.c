#include "fft_setup.h"
#ifdef USE_FFTW_WISDOM
#include "fft_wisdom.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

void setup_fftw_plans_full(int N, int narray, fftw_complex_t *plan_buffer,
                           fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *dummy_1d = NULL;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Ensure N^2 * sizeof(fftw_complex_t) is a multiple of 64 for AVX-512 alignment
    // of consecutive arrays in plan_many_dft (N divisible by 4 is sufficient)
    if (N <= 0 || (N % 4) != 0) {
        fprintf(stderr, "[ERROR] N (PPD) must be positive and divisible by 4 for FFT alignment (got N=%d)\n", N);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ====================================================================================
    // INITIALIZE FFTW THREADING 
    // ====================================================================================
    // Thread init must precede other FFTW calls (including wisdom import).

    static int fftw_threads_initialized = 0;
    if (!fftw_threads_initialized) {
        if (FFTW_INIT_THREADS() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (%s precision)\n", rank, PRECISION_NAME);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (rank == 0) {
            printf("[FFTW-THREADING] %s precision: FFTW threads init (per-plan counts: 2D=OMP_MAX, 1D=1)\n",
                   PRECISION_NAME);
        }
        fftw_threads_initialized = 1;
    }
    // ====================================================================================

#ifdef USE_FFTW_WISDOM
    // Each rank imports the same wisdom file (see wisdom_rank0).
    fft_wisdom_import_from_file(rank);

#endif
    
    // plan_buffer: one N×N complex plane (staged 2D FFT). narray is unused for 2D planning.
    (void)narray;
    if (plan_buffer == NULL) {
        fprintf(stderr, "[ERROR] Failed to provide plan_buffer for 2D FFT plan (need plane w/ PPD^2 complexes)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    {
        int fft_threads_2d = omp_get_max_threads();
        FFTW_PLAN_WITH_NTHREADS(fft_threads_2d);
        if (rank == 0) {
            printf("[FFTW-THREADING] 2D single-plane plan (staged): FFTW_PLAN_WITH_NTHREADS(%d)\n",
                   fft_threads_2d);
            fflush(stdout);
        }
    }

    /* In-place 2D complex DFT on one contiguous N×N plane (row-major). Execution loops
     * over narray planes in ZD_MPI_generation.c with copy-in/copy-out. */
    *plan_2d_out = FFTW_PLAN_DFT_2D(N, N, plan_buffer, plan_buffer, FFT_SIGN, FFTW_PLANNER_FLAGS);

    // 1D Y FFT: single FFTW thread (OpenMP parallelizes across pencils; staged buffers in z_streaming)
    FFTW_PLAN_WITH_NTHREADS(1);
    if (rank == 0) {
        printf("[FFTW-THREADING] 1D Y plan: FFTW_PLAN_WITH_NTHREADS(1)\n");
        fflush(stdout);
    }
    
    // Create 1D FFT plan (planner flags: FFTW_PLANNER_FLAGS in config.h)
    // note: ALIGN_BYTES defined in config.h to be 4096
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_PLANNER_FLAGS);
    
    free(dummy_1d);
    
    // Check plan creation
    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[ERROR] Failed to create FFT plans (one or both plans are NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

#ifdef USE_FFTW_WISDOM
    fft_wisdom_export_rank0(rank);
#endif
}

