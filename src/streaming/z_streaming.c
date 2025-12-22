// ====================================================================================
// Z-STREAMING MODULE
// ====================================================================================

#include "z_streaming.h"
#include "../types.h"  // For ZSLAB macro
#include "../config.h"  // For DEBUG_PRINTS, SKIP_VERIFICATION
#include "../precision.h"  // For real_t, fabs_t, fmax_t
#include <stdio.h>
#include <mpi.h>

// ====================================================================================
// Z-slab stream unpacking for Abcus-compatible output (?)
// Unpacks one Z-slab (cuz of mem) from the MPI receive buffer
// Applies 1D FFT along the Y-direction
// ** To-do: Put into format for zeldovich-PLT writing interface 

// ** NOTE: Batch-aware unpacking: The receive buffer is organized by batches, 
// so the function computes cumulative batch offsets to find the correct data location
// ====================================================================================

void z_streaming_unpack(
    int rank, int N, int narray,
    int z_global,                          // Which Z-slab to process
    GridBounds my_bounds,                  // My (X,Z) region bounds
    int my_pencils,                        // Total pencils in my region
    fftw_complex_t *recv_buffer,           // Source data (source-grouped)
    int64_t *recv_displs_src,              // Base offset per source
    int *src_total_slices,                 // Total slices per source (unused, kept for API)
    int *y_owner_src,                      // Y -> src mapping
    int *y_src_local_idx,                  // Y -> local_idx mapping (unused, kept for API)
    int *y_batch_idx,                      // Y -> batch mapping
    int *y_slice_idx_in_batch,             // Y -> slice_idx in batch mapping
    int **src_batch_slice_counts,          // [src][batch] -> slice count
    int global_max_batches,                // Total number of batches
    fftw_complex_t *local_z_slab,          // Destination buffer (one Z-slab)
    fftw_plan_t plan_1d_y)                 // FFT plan for Y-direction
{
    (void)rank;                // Unused
    (void)src_total_slices;    // Unused but kept for API consistency
    (void)y_src_local_idx;     // Unused but kept for API consistency
    (void)global_max_batches;  // Unused but kept for API consistency
    
    int x_count = my_bounds.x_end - my_bounds.x_start;
    int z_count = my_bounds.z_end - my_bounds.z_start;
    int z_idx = z_global - my_bounds.z_start;
    
    // Sanity check
    if (z_global < my_bounds.z_start || z_global >= my_bounds.z_end) {
        fprintf(stderr, "[ERROR] Rank %d: z_global=%d out of bounds [%d,%d)\n",
                rank, z_global, my_bounds.z_start, my_bounds.z_end);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ========== UNPACKING: Extract this Z-slab from recv_buffer ==========
    // Loop over all arrays, all X in my region, all Y
    // v14: use PERIODIC_X for actual x_idx access
    // Fixed bug: Use batch and slice_idx information to compute correct offset
    // Data arrives as [array][slice_batch_local][pencil] per batch
    // We need to compute cumulative offset of all previous batches + current batch offset
    // Write to local_z_slab in [Array][X][Y] format (Y stride-1 for FFT, better cache locality)
    #pragma omp parallel for collapse(3)
    for (int array_idx = 0; array_idx < narray; array_idx++) {
        for (int x_idx = 0; x_idx < x_count; x_idx++) {
            for (int y = 0; y < N; y++) {

                // Use pre-built lookup arrays to find which rank owns
                // this Y-slice and which batch it came from
                int src = y_owner_src[y];
                int batch = y_batch_idx[y];
                int slice_idx = y_slice_idx_in_batch[y];
                
                // Calculate pencil index in full region: (x_local, z_local)
                int pencil_idx = x_idx * z_count + z_idx;
                
                // Compute cumulative offset of all previous batches from this source
                int64_t batch_offset = 0;
                for (int b = 0; b < batch; b++) {
                    int slices_in_batch = src_batch_slice_counts[src][b];
                    batch_offset += (int64_t)slices_in_batch * my_pencils * narray;
                }
                
                // Calculate position in recv_buffer
                // Actual layout: [src][array][batch_0_slice_0][pencil], [src][array][batch_0_slice_1][pencil], ...
                // Packing order per batch: [array][slice][pencil]
                // So for this batch: offset = batch_offset + array_idx * (slices_in_this_batch * my_pencils) + slice_idx * my_pencils + pencil_idx
                int slices_in_this_batch = src_batch_slice_counts[src][batch];
                int64_t recv_offset = recv_displs_src[src]
                                    + batch_offset
                                    + (int64_t)array_idx * slices_in_this_batch * my_pencils
                                    + (int64_t)slice_idx * my_pencils
                                    + pencil_idx;
                
                // V14: Apply periodic BC when writing to local_z_slab
                // x_idx is offset from my_bounds.x_start, which can be negative
                // No periodic needed here since ZSLAB uses local x_idx (not global X)
                // The periodic wrapping was already handled during packing
                
                // Extract to local_z_slab: [Array][X][Y] format (Y stride-1 for FFT)
                ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0] = recv_buffer[recv_offset][0];
                ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1] = recv_buffer[recv_offset][1];
            }
        }
    }
    
    // ========== DEBUG: Check imaginary parts BEFORE 1D FFT ==========
    #if DEBUG_PRINTS && !SKIP_VERIFICATION
    if (rank == 0 && z_global == my_bounds.z_start) {
        real_t debug_max_real_before[4] = {0, 0, 0, 0};
        real_t debug_max_imag_before[4] = {0, 0, 0, 0};
        for (int array_idx = 0; array_idx < narray; array_idx++) {
            for (int x_idx = 0; x_idx < x_count; x_idx++) {
                for (int y = 0; y < N; y++) {
                    double re = fabs_t(ZSLAB(array_idx, x_idx, y, N, narray, x_count)[0]);
                    double im = fabs_t(ZSLAB(array_idx, x_idx, y, N, narray, x_count)[1]);
                    debug_max_real_before[array_idx] = fmax_t(debug_max_real_before[array_idx], re);
                    debug_max_imag_before[array_idx] = fmax_t(debug_max_imag_before[array_idx], im);
                }
            }
        }
        printf("[DEBUG Z=%d] Before 1D FFT (after unpack) - Max imag per array: ", z_global);
        for (int a = 0; a < narray; a++) {
            printf("A%d=%.3e ", a, debug_max_imag_before[a]);
        }
        printf("\n");
    }
    #endif
    
    // ========== 1D FFT: Apply along Y-direction for each (Array, X) ==========
    // In [Array][X][Y] format, Y is stride-1 for fixed (array_idx, x_idx)
    #pragma omp parallel for collapse(2)
    for (int array_idx = 0; array_idx < narray; array_idx++) {
        for (int x_idx = 0; x_idx < x_count; x_idx++) {
            fftw_complex_t *y_data = &ZSLAB(array_idx, x_idx, 0, N, narray, x_count);
            FFTW_EXECUTE_DFT(plan_1d_y, y_data, y_data);
        }
    }
    
    // local_z_slab now contains FFT'd data for this Z-slab, ready for writing
}

