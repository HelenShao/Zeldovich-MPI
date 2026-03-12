#include "fft_setup.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

// todo: add FFTW_WISDOM?
// fftw_import_wisdom_file("wisdom_double.txt")
// Export wisdom after creating plans

void setup_fftw_plans_full(int N, int narray, fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *dummy_2d = NULL;
    fftw_complex_t *dummy_1d = NULL;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // ====================================================================================
    // INITIALIZE FFTW THREADING 
    // ====================================================================================
    // FFTW threading must be initialized before creating plans (to use OMP)

    static int fftw_threads_initialized = 0;
    if (!fftw_threads_initialized) {
        int nthreads = omp_get_max_threads();
        
        // Use all OpenMP threads inside FFTW for each plan execution.
        // We no longer split threads between outer narray parallelism and inner FFTW threading.
        int fft_threads = nthreads;
        
        // Using macros from precision.h for double or single 
        if (FFTW_INIT_THREADS() == 0) {
            fprintf(stderr, "[ERROR] Rank %d: Failed to initialize FFTW threads (%s precision)\n", rank, PRECISION_NAME);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        FFTW_PLAN_WITH_NTHREADS(fft_threads);
        if (rank == 0) {
            printf("[FFTW-THREADING] %s precision: Using %d FFTW threads per FFT (narray=%d, total_threads=%d)\n", 
                   PRECISION_NAME, fft_threads, narray, nthreads);
        }
        
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
    
    *plan_2d_out = FFTW_PLAN_DFT_2D(N, N, dummy_2d, dummy_2d, 
                                     FFT_SIGN, FFTW_MEASURE);
    
    free(dummy_2d);
    
    // Create 1D FFT plan with FFTW_MEASURE
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_MEASURE);
    
    free(dummy_1d);
    
    // Ceck plan creation
    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[ERROR] Failed to create FFT plans (one or both plans are NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

