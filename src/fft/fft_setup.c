// ====================================================================================
// HERMITIAN 3D MATRIX MPI - FFT SETUP MODULE
// ====================================================================================

#include "fft_setup.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

// ====================================================================================
// FUNCTION IMPLEMENTATIONS
// ====================================================================================

void setup_fftw_plans_full(int N, fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out)
{
    fftw_complex_t *dummy_2d = NULL;
    fftw_complex_t *dummy_1d = NULL;
    
    // Create 2D FFT plan (for Y-slice processing)
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
    
    // Create 1D FFT plan (for Y-direction FFT on pencils, used in Stage 6)
    if (posix_memalign((void**)&dummy_1d, ALIGN_BYTES, 
                       sizeof(fftw_complex_t) * N) != 0) {
        fprintf(stderr, "[ERROR] Failed to allocate dummy_1d for 1D FFT plan creation\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    *plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, 
                                     FFT_SIGN, FFTW_ESTIMATE);
    
    // CRITICAL: Free dummy memory immediately after plan creation
    free(dummy_1d);
    
    // Verify plan creation succeeded
    if (*plan_2d_out == NULL || *plan_1d_out == NULL) {
        fprintf(stderr, "[ERROR] Failed to create FFT plans (one or both plans are NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

