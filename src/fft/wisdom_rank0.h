#ifndef WISDOM_RANK0_H
#define WISDOM_RANK0_H

#include "../precision.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Import/export wisdom at FFTW_WISDOM_FILENAME (see fft_wisdom.h),
 * initialize FFTW threads, create 2D DFT_2D + 1D DFT plans on plan_buffer (in-place 2D on one N×N plane),
 * export accumulated wisdom to the same path.
 *
 * plan_buffer must be posix_memalign(ALIGN_BYTES, N * N * sizeof(fftw_complex_t)).
 * narray is accepted for CLI compatibility; 2D wisdom does not depend on it.
 * Caller destroys returned plans after this returns.
 */
int wisdom_rank0_plans_and_export(int N, int narray, fftw_complex_t *plan_buffer,
                                  fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out);

#ifdef __cplusplus
}
#endif

#endif
