#ifndef KERNEL_PASS_H
#define KERNEL_PASS_H

// ====================================================================================
// ONE KERNEL TERM:  out += coef * J(k) * Q(k),  Q = FFT[ IFFT(G*phi) * IFFT(H*phi) ]
// ====================================================================================
//   Stage A (forward exchange, batched over y pairs):
//       z(k) = (G + iH) * phi at padded positions -> 2D bwd FFT per y-slice
//       -> pack_slices_to_send_buffer -> Alltoallv into the persistent recv_buffer
//   Stage B (per z-slab of my pencils):
//       unpack (zero for padded y rows) -> 1D bwd FFT in y -> Q = Im(z^2)/2/M^3
//       -> 1D fwd FFT in y -> keep only y_p in [0, N/2] (y crop) in q_half
//       (q_half lives in recv_buffer: see kernel_pass_init)
//   Stage C (return exchange, batched over y pairs):
//       q_half column y_p -> Alltoallv to the y_p owner -> M x M plane -> 2D fwd FFT
//       -> read the N x N original modes (x,z crop) -> out += coef * J * Q
// ====================================================================================

#include <mpi.h>
#include "precision.h"
#include "kernel_grid.h"
#include "kernel_terms.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const KernelLayout *K;
    MPI_Comm comm;
    double kfund;                 // 2*pi / BoxSize

    fftw_complex_t **phi;         // [my_num_pairs] primary slices, N*N [z][x]
    fftw_complex_t **out;         // [my_num_pairs] accumulated output, N*N [z][x]

    fftw_complex_t *recv_buffer;  // K->recv_total = N * my_pencils, [slot][pencil]
    fftw_complex_t *q_half;       // [(N/2+1)][my_pencils], aliases recv_buffer (not owned)

    fftw_plan_t plan2d_bwd, plan2d_fwd;  // M x M in-place
    fftw_plan_t plan1d_bwd, plan1d_fwd;  // length M in-place

    // Accumulated wall times (s)
    double t_fill_fft2d, t_pack, t_comm_fwd;
    double t_zslab;
    double t_comm_ret, t_fft2d_accum;
} KernelPass;

void kernel_pass_init(KernelPass *P, const KernelLayout *K, MPI_Comm comm, double box_size,
                      fftw_complex_t **phi, fftw_complex_t **out);
void kernel_pass_free(KernelPass *P);

void run_kernel_term(KernelPass *P, const KernelShape *shape, double coef);

// Aligned allocation that aborts on failure (size rounded up to a page)
fftw_complex_t *kernel_alloc_complex(int64_t count, const char *what);

#ifdef __cplusplus
}
#endif

#endif
