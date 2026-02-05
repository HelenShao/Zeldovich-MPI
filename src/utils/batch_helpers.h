#ifndef HERMITIAN_BATCH_HELPERS_H
#define HERMITIAN_BATCH_HELPERS_H

#include "config.h"
#include "types.h"
#include <mpi.h>

// ====================================================================================
// Calculate how many Y-slices a given rank processes in a given batch
// Returns: 0 (no slices), 1 (self-conjugate), or 2 (conjugate pair)
// This allows each rank to independently calculate send/recv counts per batch

int get_rank_batch_slice_count(int target_rank, int batch_idx, int N, int num_ranks);

// ====================================================================================

// Get the specific Y-values that a rank is processing in a given batch
// Fills out_y_values array (must be pre-allocated with size >= 2)
// Sets out_count to number of Y-values (0, 1, or 2)

void get_rank_batch_y_values(int target_rank, int batch_idx, int N, int num_ranks,
                              int *out_y_values, int *out_count);
// ====================================================================================

// Calculate sendcounts, recvcounts,  total send and recv counts, and displacements for a given batch
// Allocates new arrays (caller must free them)

void calculate_batch_send_recv_counts(
    int rank, int num_ranks, int N, int narray, int batch_idx,
    int my_batch_slice_count, int my_pencils,
    int **out_sendcounts, int **out_sdispls,
    int **out_recvcounts, int **out_rdispls,
    int *out_total_send, int *out_total_recv);
// ====================================================================================

#endif 

