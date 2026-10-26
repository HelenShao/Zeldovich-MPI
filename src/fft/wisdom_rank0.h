#ifndef WISDOM_RANK0_H
#define WISDOM_RANK0_H

#include "../precision.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Create 2D + 1D FFTW plans on plan_buffer.
 * plan_buffer must be posix_memalign(ALIGN_BYTES, N * N * sizeof(fftw_complex_t)).
 * @param save_to_file Non-zero to also write wisdom to zd_wisdom_rank0_file(); zero to keep wisdom in memory only.
 */
int wisdom_rank0_plans(int N, int narray, fftw_complex_t *plan_buffer, fftw_plan_t *plan_2d_out,
                       fftw_plan_t *plan_1d_out, int save_to_file);

#ifdef __cplusplus
}
#endif

#endif
