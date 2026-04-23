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
    int fft_threads_2d = 1;
    if (!fftw_threads_initialized) {
        int nthreads = omp_get_max_threads();
        
        // All threads for FFTW (no outer omp over narray).
        fft_threads_2d = nthreads;
        
        if (FFTW_INIT_THREADS() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (%s precision)\n", rank, PRECISION_NAME);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        FFTW_PLAN_WITH_NTHREADS(fft_threads_2d);
        if (rank == 0) {
            printf("[FFTW-THREADING] %s precision: 2D plan uses %d FFTW threads; 1D plan will use 1 thread\n",
                   PRECISION_NAME, fft_threads_2d);
        }
        
        fftw_threads_initialized = 1;
    } else {
        fft_threads_2d = omp_get_max_threads();
    }
    // ====================================================================================

#ifdef USE_FFTW_WISDOM
    // Each rank imports the same wisdom file (see wisdom_rank0).
    fft_wisdom_import_from_file(rank);

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

    // Verify thread count for 2D planning/execution.
    if (rank == 0) {
        int planner_n = omp_get_max_threads();
        printf("[FFTW-THREADING] 2D planner thread target (OMP max): %d\n",
               planner_n);
        fflush(stdout);
    }
    
    // Force 1 thread for 1D plan creation (and later execution of this plan handle).
    FFTW_PLAN_WITH_NTHREADS(1);

    // Create 1D FFT plan with FFTW_MEASURE
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_MEASURE);

    if (rank == 0) {
        printf("[FFTW-THREADING] 1D plan created with FFTW threads: 1\n");
        fflush(stdout);
    }
    
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

