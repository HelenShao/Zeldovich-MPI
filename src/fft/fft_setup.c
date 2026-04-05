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
    fftw_complex_t *buf_2d = NULL;
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
        int nthreads = omp_get_max_threads();
        
        // All threads for FFTW (no outer omp over narray).
        int fft_threads = nthreads;
        
        if (FFTW_INIT_THREADS() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (%s precision)\n", rank, PRECISION_NAME);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        FFTW_PLAN_WITH_NTHREADS(fft_threads);
        if (rank == 0) {
            printf("[FFTW-THREADING] %s precision: %d FFTW threads (no outer OMP over narray)\n", 
                   PRECISION_NAME, fft_threads);
        }
        
        fftw_threads_initialized = 1;
    }
    // ====================================================================================

#ifdef USE_FFTW_WISDOM
    // Load per-rank wisdom first so all FFTW_MEASURE calls below (probe on rank 0,
    // then real 2D/1D plans) can consult it.
    fft_wisdom_import_per_rank(rank, MPI_COMM_WORLD);

    // Rank 0: optional sanity check that MEASURE still adds bytes vs current store.
    // Do not FFTW_FORGET_WISDOM here — that would drop imported wisdom before planning.
    if (rank == 0) {
        char *pre = FFTW_EXPORT_WISDOM_TO_STRING();
        size_t pre_len = pre ? strlen(pre) : 0;
        if (pre) FFTW_FREE(pre);

        fftw_complex_t *probe = NULL;
        if (posix_memalign((void**)&probe, 64, sizeof(fftw_complex_t) * 16) == 0) {
            fftw_plan_t p = FFTW_PLAN_DFT_1D(16, probe, probe, FFT_SIGN, FFTW_MEASURE);
            if (p) FFTW_DESTROY_PLAN(p);
            free(probe);
        }

        char *post = FFTW_EXPORT_WISDOM_TO_STRING();
        size_t post_len = post ? strlen(post) : 0;
        if (post) FFTW_FREE(post);

        printf("[FFTW-WISDOM-DIAG] Probe (after import): wisdom before=%zu after=%zu bytes "
               "(delta=%zd). %s\n",
               pre_len, post_len, (ssize_t)(post_len - pre_len),
               post_len > pre_len ? "OK - MEASURE appended to wisdom store."
                                  : "No growth (often OK if probe plan already in imported wisdom).");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    // plan_buffer must be provided before setup
    if (plan_buffer == NULL) {
        fprintf(stderr, "[ERROR] Failed to provide plan_buffer for 2D batched FFT plan\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    buf_2d = plan_buffer;
    
    /* Memory layout for plan_many_dft below: narray contiguous N×N complex planes.
     * Plane index a (0 <= a < narray) starts at buf_2d + a * N * N.
     * Within a plane, row-major C order: element (i,j) at buf_2d[a*N*N + i*N + j]. */
    {
        int n[2] = { N, N };
        *plan_2d_out = FFTW_PLAN_MANY_DFT(
            2,      // 2D transform
            n,      // each transform is size N×N
            narray, // number of transforms in the batch
            buf_2d, NULL, 1, N * N, // input data
            buf_2d, NULL, 1, N * N, // output data
            FFT_SIGN, // FFT direction
            FFTW_MEASURE);
    }

    // Verify thread count that FFTW will use when executing this plan (must match plan_with_nthreads)
    if (rank == 0) {
        int planner_n = FFTW_PLANNER_NTHREADS();
        printf("[FFTW-THREADING] Planner nthreads (used at execute): %d\n", planner_n);
        fflush(stdout);
    }
    
    // Create 1D FFT plan with FFTW_MEASURE
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_MEASURE);
    
    free(dummy_1d);
    
    // Check plan creation
    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[ERROR] Failed to create FFT plans (one or both plans are NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

#ifdef USE_FFTW_WISDOM
    // each rank exports its own wisdom file for reuse.
    fft_wisdom_export_per_rank(rank);
#endif
}

