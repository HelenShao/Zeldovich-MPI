// ====================================================================================
// FFT SETUP MODULE
// ====================================================================================

#include "fft_setup.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

void setup_fftw_plans_full(int N, fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *dummy_2d = NULL;
    fftw_complex_t *dummy_1d = NULL;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // ====================================================================================
    // INITIALIZE FFTW THREADING (CRITICAL FOR PERFORMANCE)
    // ====================================================================================
    // FFTW threading must be initialized before creating plans.
    // This enables FFTW to use OpenMP threads for FFT computation.
    // Without this, FFTW only uses ONE thread regardless of OMP_NUM_THREADS setting!
    static int fftw_threads_initialized = 0;
    if (!fftw_threads_initialized) {
        int nthreads = omp_get_max_threads();
        
        #ifdef USE_DOUBLE_PRECISION
        // Double precision
        if (fftw_init_threads() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (double precision)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fftw_plan_with_nthreads(nthreads);
        if (rank == 0) {
            printf("[FFTW-THREADING] Double precision: Initialized with %d threads per rank\n", nthreads);
        }
        #else
        // Single precision (default)
        if (fftwf_init_threads() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (single precision)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fftwf_plan_with_nthreads(nthreads);
        if (rank == 0) {
            printf("[FFTW-THREADING] Single precision: Initialized with %d threads per rank\n", nthreads);
        }
        #endif
        
        fftw_threads_initialized = 1;
    }
    // ====================================================================================
    
    // Create 2D FFT plan (for XY-direction within Y-slices)
    if (posix_memalign((void**)&dummy_2d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_2d for 2D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_2d_out = FFTW_PLAN_DFT_2D(N, N, dummy_2d, dummy_2d, 
                                     FFT_SIGN, FFTW_ESTIMATE);
    
    // CRITICAL: Free dummy memory immediately after plan creation
    // The plan is now independent of the dummy buffer
    free(dummy_2d);
    
    // Create 1D FFT plan (for Y-direction FFT on pencils)
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_ESTIMATE);
    
    // CRITICAL: Free dummy memory immediately after plan creation
    free(dummy_1d);
    
    // Verify plan creation
    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[ERROR] Failed to create FFT plans (one or both plans are NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

