#ifndef FFT_SETUP_H
#define FFT_SETUP_H

#include "../config.h"
#include "../precision.h"
#include "../types.h"

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// Create FFTW plans for 2D and 1D transforms
// - plan_buffer: Buffer for 2D plan: N * N fftw_complex_t (one X–Z plane). Aligned (e.g. posix_memalign ALIGN_BYTES).
// - plan_2d_out: Single-plane in-place FFTW_PLAN_DFT_2D(N,N,...) on plan_buffer.
//   Planned with FFTW_PLAN_WITH_NTHREADS(omp_get_max_threads()). narray is passed for API compatibility only.
// - plan_1d_out: Y-direction 1D FFT (size N). Planned with FFTW_PLAN_WITH_NTHREADS(1);
//   execution uses staged per-OMP-thread buffers in z_streaming (no nested FFTW threads).
// ====================================================================================

void setup_fftw_plans_full(int N, int narray, fftw_complex_t *plan_buffer,
                           fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out);

#ifdef __cplusplus
}
#endif

#endif 

