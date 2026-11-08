// ====================================================================================
// KERNEL GRID LAYOUT
// ====================================================================================

#include "kernel_grid.h"
#include "mpi_topology.h"
#include "utils/batch_helpers.h"
#include "utils/decomposition.h"
#include <stdio.h>
#include <stdlib.h>

int kernel_choose_pad_M(int N, int compat_2N)
{
    if (compat_2N) return 2 * N;
    int M = (3 * N + 1) / 2;
    while ((M % 4) != 0 || ((M - N) % 2) != 0) M++;
    return M;
}

GridBounds kernel_pencil_bounds(int dest, int M, int num_ranks, int size_x, int size_z)
{
    // cpd = M makes the CPD-aligned split a plain integer split of [0, M). The same
    // call is made inside pack_slices_to_send_buffer(), so send and recv agree.
    ExtendedGridBounds ext = get_extended_grid_bounds_CPD_aligned(dest, M, num_ranks, size_x, size_z, M);
    return ext.core;
}

void kernel_layout_init(KernelLayout *K, int N, int M, int num_ranks, int rank,
                        int size_x, int size_z)
{
    K->N = N;
    K->M = M;
    K->num_ranks = num_ranks;
    K->rank = rank;
    K->size_x = size_x;
    K->size_z = size_z;

    // ---------- Pencil decomposition ----------
    K->bounds = (GridBounds*)malloc(sizeof(GridBounds) * num_ranks);
    for (int d = 0; d < num_ranks; d++) {
        K->bounds[d] = kernel_pencil_bounds(d, M, num_ranks, size_x, size_z);
    }
    K->my_bounds = K->bounds[rank];
    K->x_count = K->my_bounds.x_end - K->my_bounds.x_start;
    K->z_count = K->my_bounds.z_end - K->my_bounds.z_start;
    K->my_pencils = kernel_region_size(&K->my_bounds);

    // ---------- Pair ownership (same distribution as batch_helpers.c) ----------
    int total_pairs = N / 2 + 1;
    K->num_batches = (total_pairs + num_ranks - 1) / num_ranks;
    K->batch_y_p = (int*)malloc(sizeof(int) * (size_t)num_ranks * K->num_batches);
    K->src_batch_slice_counts = (int**)malloc(sizeof(int*) * num_ranks);
    K->my_num_pairs = 0;

    K->y_base = (int64_t*)malloc(sizeof(int64_t) * M);
    int *y_seen = (int*)calloc(N, sizeof(int));
    for (int y = 0; y < M; y++) K->y_base[y] = -1;

    // ---------- Forward recv layout: [src][batch][slice][pencil] ----------
    K->recv_displs_src = (int64_t*)malloc(sizeof(int64_t) * num_ranks);
    int64_t offset = 0;
    for (int src = 0; src < num_ranks; src++) {
        K->recv_displs_src[src] = offset;
        K->src_batch_slice_counts[src] = (int*)calloc((size_t)K->num_batches, sizeof(int));
        for (int b = 0; b < K->num_batches; b++) {
            int y_values[2];
            int count = 0;
            get_rank_batch_y_values(src, b, N, num_ranks, y_values, &count);
            K->src_batch_slice_counts[src][b] = count;
            K->batch_y_p[(int64_t)src * K->num_batches + b] = (count > 0) ? y_values[0] : -1;
            if (src == rank && count > 0) K->my_num_pairs++;

            for (int s = 0; s < count; s++) {
                int y = y_values[s];
                if (y < 0 || y >= N || y_seen[y]) {
                    fprintf(stderr, "[Kernel_MPI] ERROR: y=%d assigned twice or out of range (src %d, batch %d)\n",
                            y, src, b);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                y_seen[y] = 1;
                K->y_base[pad_idx(y, N, M)] = offset + (int64_t)s * K->my_pencils;
            }
            offset += (int64_t)count * K->my_pencils;
        }
    }
    K->recv_total = offset;

    for (int y = 0; y < N; y++) {
        if (!y_seen[y]) {
            fprintf(stderr, "[Kernel_MPI] ERROR: original y=%d not owned by any rank\n", y);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    free(y_seen);
}

void kernel_layout_free(KernelLayout *K)
{
    for (int src = 0; src < K->num_ranks; src++) free(K->src_batch_slice_counts[src]);
    free(K->src_batch_slice_counts);
    free(K->batch_y_p);
    free(K->bounds);
    free(K->recv_displs_src);
    free(K->y_base);
}
