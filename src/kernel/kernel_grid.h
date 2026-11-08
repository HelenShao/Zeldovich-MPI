#ifndef KERNEL_GRID_H
#define KERNEL_GRID_H

// ====================================================================================
// KERNEL GRID LAYOUT
// ====================================================================================
// Two grids are in play:
//   - the original N^3 grid (files, y-pair ownership, G/H/J evaluation, output)
//   - the padded M^3 grid (FFTs, pencil decomposition)
//
// Padding and cropping are never done as separate steps: original index i maps to
// padded index pad_idx(i), the same mapping as fftshift -> np.pad -> ifftshift in
// kernel_MPI.py (frequency n mod N -> n mod M). Loops over the original N indices
// write (pad) or read (crop) at pad_idx(); the rest of the padded grid is untouched.
//
// Y-slice ownership (which rank reads/writes which y_primary) is the file ownership
// from batch_helpers.c on the ORIGINAL grid: N/2+1 conjugate pairs split over ranks.
// Batch b on rank r processes r's b-th pair (y_p, and y_m = N - y_p if not self-conj).
//
// Pencil ownership is the (x,z) block of the PADDED grid from decomposition.c
// (get_extended_grid_bounds_CPD_aligned with cpd = M -> plain integer split).
// ====================================================================================

#include <stdint.h>
#include "config.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

// Signed frequency of FFT index i on an n-point grid (np.fft.fftfreq(n) * n)
static inline int kfreq(int i, int n) { return i < n / 2 ? i : i - n; }

// Original-grid index -> padded-grid index (same frequency)
static inline int pad_idx(int i, int N, int M) { return i < N / 2 ? i : i + (M - N); }

// Smallest valid padded size >= target: M % 4 == 0 (fft_setup.c) and (M - N) even
int kernel_choose_pad_M(int N, int compat_2N);

// Pencil (x,z) bounds of rank dest on the padded grid
GridBounds kernel_pencil_bounds(int dest, int M, int num_ranks, int size_x, int size_z);

typedef struct {
    int N, M, num_ranks, rank, size_x, size_z;

    // Pencil decomposition (padded grid)
    GridBounds *bounds;          // [num_ranks]
    GridBounds my_bounds;
    int x_count, z_count;
    int64_t my_pencils;

    // Pair ownership (original grid)
    int num_batches;             // max pairs per rank (= global batch count)
    int my_num_pairs;
    int *batch_y_p;              // [num_ranks * num_batches], y_primary or -1
    int **src_batch_slice_counts;// [num_ranks][num_batches], 0/1/2

    // Forward recv_buffer layout: [src][batch][slice][pencil] (narray = 1)
    int64_t *recv_displs_src;    // [num_ranks]
    int64_t recv_total;

    // Padded-y lookup for unpacking: y_base[y] = recv_buffer offset of pencil 0
    // for padded row y, or -1 for a zero (padding) row
    int64_t *y_base;             // [M]
} KernelLayout;

void kernel_layout_init(KernelLayout *K, int N, int M, int num_ranks, int rank,
                        int size_x, int size_z);
void kernel_layout_free(KernelLayout *K);

// y_primary that rank r processes in batch b, or -1
static inline int kernel_batch_y_p(const KernelLayout *K, int r, int b)
{
    return K->batch_y_p[(int64_t)r * K->num_batches + b];
}

static inline int64_t kernel_region_size(const GridBounds *g)
{
    return (int64_t)(g->x_end - g->x_start) * (int64_t)(g->z_end - g->z_start);
}

#ifdef __cplusplus
}
#endif

#endif
