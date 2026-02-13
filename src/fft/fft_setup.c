#include "fft_setup.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

// todo: add FFTW_WISDOM
// Double: fftw_import_wisdom_file("wisdom_double.txt")
// Export wisdom after creating plans

void setup_fftw_plans_full(int N, fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *dummy_2d = NULL;
    fftw_complex_t *dummy_1d = NULL;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // ====================================================================================
    // INITIALIZE FFTW THREADING
    // ====================================================================================
    // All plans are created with nthreads=1 (single-threaded) so that they are SAFE
    // for concurrent execution from multiple OMP threads via FFTW_EXECUTE_DFT.
    //
    // Why single-threaded plans + outer OMP parallelism?
    //   1. FFTW plans with nthreads>1 have SHARED internal scratch buffers.
    //      Calling FFTW_EXECUTE_DFT concurrently from multiple threads with the
    //      same multi-threaded plan causes DATA RACES on those scratch buffers.
    //   2. Outer OMP parallelism (many independent single-threaded FFTs) outperforms
    //      inner FFTW threading when narray * x_count >> num_threads (always here).
    //   3. Zero sync overhead between FFTs; better cache locality per thread.
    // ====================================================================================
    static int fftw_threads_initialized = 0;
    if (!fftw_threads_initialized) {
        int nthreads = omp_get_max_threads();
        
        #ifdef USE_DOUBLE_PRECISION
        if (fftw_init_threads() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (double precision)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fftw_plan_with_nthreads(1);
        if (rank == 0) {
            printf("[FFTW-THREADING] Double precision: Plans use 1 thread (OMP_NUM_THREADS=%d for outer parallelism)\n", nthreads);
        }
        #else
        if (fftwf_init_threads() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (single precision)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fftwf_plan_with_nthreads(1);
        if (rank == 0) {
            printf("[FFTW-THREADING] Single precision: Plans use 1 thread (OMP_NUM_THREADS=%d for outer parallelism)\n", nthreads);
        }
        #endif
        
        fftw_threads_initialized = 1;
    }
    // ====================================================================================
    
    // Create 2D FFT plan with FFTW_MEASURE for better algorithm selection.
    // Extra planning time is amortized over thousands of executions.
    if (posix_memalign((void**)&dummy_2d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_2d for 2D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0) {
        printf("[FFTW-PLAN] Creating 2D plan (%d x %d) with FFTW_MEASURE...\n", N, N);
    }
    *plan_2d_out = FFTW_PLAN_DFT_2D(N, N, dummy_2d, dummy_2d, 
                                     FFT_SIGN, FFTW_MEASURE);
    
    free(dummy_2d);
    
    // Create 1D FFT plan with FFTW_MEASURE
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0) {
        printf("[FFTW-PLAN] Creating 1D plan (N=%d) with FFTW_MEASURE...\n", N);
    }
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_MEASURE);
    
    free(dummy_1d);
    
    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[ERROR] Failed to create FFT plans (one or both plans are NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0) {
        printf("[FFTW-PLAN] All plans created (single-threaded, FFTW_MEASURE)\n");
    }
}

