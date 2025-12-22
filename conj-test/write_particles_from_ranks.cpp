// ====================================================================================
// WRITE PARTICLES FROM PER-RANK Z-SLABS (NO REASSEMBLY)
// ====================================================================================
// Reads per-rank i-slab files written by main.cpp and calls WriteParticlesSlab_range
// to write particle data directly (without reassembly).
//
// This program:
// 1. Reads command-line arguments
// 2. Loads or creates simulation parameters
// 3. Initializes output buffers
// 4. For each rank and i-slab:
//    - Reads the per-rank binary file: rank_{rank}/i{i}_slab_N{N}.bin
//    - Determines rank's X-range using grid decomposition
//    - Extracts pointers to the 4 arrays
//    - Calls WriteParticlesSlab_range to write per-rank particle data
//
// Example:
//   ./write_particles_from_ranks . 16 4 param_N16.par 0 16
// ====================================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <vector>
#include <string>
#include <filesystem>
#include <iostream>

// Include zeldovich-PLT headers
#include <output.h>
#include <parameters.h>
#include <omp.h>

// Include grid decomposition utilities
extern "C" {
#include "utils/decomposition.h"
}

namespace fs = std::filesystem;
using Complx = std::complex<double>;

// ====================================================================================
// READ PER-RANK I-SLAB
// ====================================================================================

std::vector<Complx> read_rank_i_slab(
    const std::string& output_dir,
    int rank,
    int i,
    int N,
    int narray,
    int x_count_expected
) {
    std::vector<Complx> local_slab;
    
    // Build filename
    std::string filename = output_dir + "/rank_" + std::to_string(rank) + 
                          "/i" + std::to_string(i) + "_slab_N" + std::to_string(N) + ".bin";
    
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open file: %s\n", filename.c_str());
        return local_slab;
    }
    
    // Get file size
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    
    // Expected size: narray * N * x_count * 16 bytes (double precision complex)
    size_t expected_size = narray * N * x_count_expected * 16;
    if (file_size != (long)expected_size) {
        fprintf(stderr, "WARNING: File size mismatch for %s: expected %zu, got %ld\n",
                filename.c_str(), expected_size, file_size);
        // Try to read anyway
    }
    
    // Infer x_count from file size
    int x_count = file_size / (narray * N * 16);
    
    // Allocate and read
    local_slab.resize(narray * N * x_count);
    size_t read_size = fread(local_slab.data(), 16, narray * N * x_count, fp);
    fclose(fp);
    
    if (read_size != narray * N * x_count) {
        fprintf(stderr, "ERROR: Read %zu elements, expected %d\n", 
                read_size, narray * N * x_count);
        local_slab.clear();
    }
    
    return local_slab;
}

// ====================================================================================
// MAIN PROGRAM
// ====================================================================================

int main(int argc, char* argv[]) {
    // Force single OpenMP thread to avoid segfault in zeldovich-PLT initialization
    omp_set_num_threads(1);
    
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <output_dir> <N> <num_ranks> [param_file] [i_start] [i_end]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  output_dir: Directory containing rank_*/ subdirectories\n");
        fprintf(stderr, "  N: Resolution\n");
        fprintf(stderr, "  num_ranks: Number of MPI ranks\n");
        fprintf(stderr, "  param_file: Optional parameter file (default: create minimal)\n");
        fprintf(stderr, "  i_start: First i-slab to process (default: 0)\n");
        fprintf(stderr, "  i_end: Last i-slab to process (default: N)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Example:\n");
        fprintf(stderr, "  %s . 16 4 param_N16.par 0 16\n", argv[0]);
        return 1;
    }
    
    std::string output_dir = argv[1];
    int N = std::stoi(argv[2]);
    int num_ranks = std::stoi(argv[3]);
    int narray = 4;
    
    std::string param_file = "";
    int i_start = 0;
    int i_end = N;
    
    if (argc >= 5) {
        param_file = argv[4];
    }
    if (argc >= 6) {
        i_start = std::stoi(argv[5]);
    }
    if (argc >= 7) {
        i_end = std::stoi(argv[6]);
    }
    
    printf("================================================================================\n");
    printf("WRITE PARTICLES FROM PER-RANK Z-SLABS (NO REASSEMBLY)\n");
    printf("================================================================================\n");
    printf("Output directory: %s\n", output_dir.c_str());
    printf("Grid size N: %d\n", N);
    printf("Number of MPI ranks: %d\n", num_ranks);
    printf("Number of arrays: %d\n", narray);
    printf("i-slab range: [%d, %d)\n", i_start, i_end);
    printf("Parameter file: %s\n", param_file.empty() ? "(create minimal)" : param_file.c_str());
    printf("================================================================================\n\n");
    
    // Create Parameters object
    Parameters* param = nullptr;
    if (!param_file.empty() && fs::exists(param_file)) {
        try {
            param = new Parameters(fs::path(param_file));
            printf("Loaded parameters from: %s\n", param_file.c_str());
            printf("  ppd: %ld\n", param->ppd);
            printf("  boxsize: %f\n", param->boxsize);
            printf("  separation: %f\n", param->separation);
        } catch (const std::exception& e) {
            fprintf(stderr, "ERROR: Failed to load parameter file: %s\n", e.what());
            fprintf(stderr, "Creating minimal parameters...\n");
            param = nullptr;
        }
    }
    
    // Create minimal Parameters if needed
    if (!param) {
        std::string tmp_param_file = "/tmp/write_particles_ranks_params.par";
        FILE* tmp_fp = fopen(tmp_param_file.c_str(), "w");
        if (tmp_fp) {
            fprintf(tmp_fp, "BoxSize = 1000.0\n");
            fprintf(tmp_fp, "ZD_Version = 2\n");
            fprintf(tmp_fp, "ZD_PPD = %d\n", N);
            fprintf(tmp_fp, "ZD_Pk_scale = 1.0\n");
            fprintf(tmp_fp, "NP = %d\n", N * N * N);
            fprintf(tmp_fp, "ZD_NumBlock = 2\n");
            fprintf(tmp_fp, "CPD = %d\n", N);
            fprintf(tmp_fp, "ZD_Pk_norm = 1.0\n");
            fprintf(tmp_fp, "ZD_Pk_sigma = 1.0\n");
            fprintf(tmp_fp, "ZD_Pk_smooth = 0.0\n");
            fprintf(tmp_fp, "ZD_Pk_powerlaw_index = -2.0\n");
            fprintf(tmp_fp, "InitialConditionsDirectory = \"./output_particles\"\n");
            fprintf(tmp_fp, "InitialRedshift = 0.0\n");
            fprintf(tmp_fp, "ICFormat = \"RVZel\"\n");
            fprintf(tmp_fp, "ZD_Seed = 12345\n");
            fprintf(tmp_fp, "ZD_qPLT = 0\n");
            fprintf(tmp_fp, "ZD_qdensity = 1\n");
            fprintf(tmp_fp, "ZD_qascii = 0\n");
            fclose(tmp_fp);
            
            try {
                param = new Parameters(fs::path(tmp_param_file));
                printf("Created minimal parameters (ppd=%ld, boxsize=%f)\n", param->ppd, param->boxsize);
            } catch (const std::exception& e) {
                fprintf(stderr, "ERROR: Failed to create minimal parameters: %s\n", e.what());
                return 1;
            }
        } else {
            fprintf(stderr, "ERROR: Cannot create temporary parameter file\n");
            return 1;
        }
    }
    
    // Initialize output buffers (required by WriteParticlesSlab_range)
    printf("Initializing output...\n");
    SetupOutputDir(*param);
    double buffer_size = InitOutputBuffers(*param);
    printf("  Output directory: %s\n", param->output_dir.c_str());
    printf("  Buffer size: %.3f GiB\n", buffer_size);
    
    printf("\nProcessing per-rank i-slabs...\n");
    
    // Determine grid factors
    int grid_x, grid_z;
    calculate_grid_factors(num_ranks, &grid_x, &grid_z);
    printf("Grid decomposition: %d x %d = %d ranks\n", grid_x, grid_z, num_ranks);
    
    // Process each i-slab and rank
    int processed = 0;
    for (int i = i_start; i < i_end; i++) {
        printf("Processing i=%d...\n", i);
        
        for (int rank = 0; rank < num_ranks; rank++) {
            // Get X-range for this rank
            GridBounds bounds = get_grid_bounds(rank, N, num_ranks);
            int k_start_global = bounds.x_start;
            int k_extent = bounds.x_end - bounds.x_start;
            
            // Read per-rank i-slab
            std::vector<Complx> local_slab = read_rank_i_slab(
                output_dir, rank, i, N, narray, k_extent
            );
            
            if (local_slab.empty()) {
                printf("  Rank %d, i=%d: Skipping (file not found or error)\n", rank, i);
                continue;
            }
            
            // Verify we have the right amount of data
            if (local_slab.size() != (size_t)(narray * N * k_extent)) {
                fprintf(stderr, "ERROR: Rank %d, i=%d wrong size: %zu != %d\n",
                        rank, i, local_slab.size(), narray * N * k_extent);
                continue;
            }
            
            // Create 2D slab pointers for WriteParticlesSlab_range
            // Layout: [narray][Y][x_local]
            Complx* slab1 = &local_slab[0 * N * k_extent];  // Array 0: D + i*F
            Complx* slab2 = &local_slab[1 * N * k_extent];  // Array 1: G + i*H
            Complx* slab3 = &local_slab[2 * N * k_extent];  // Array 2: X-velocity
            Complx* slab4 = &local_slab[3 * N * k_extent];  // Array 3: Y-velocity + Z-velocity
            
            // Call WriteParticlesSlab_range
            try {
                WriteParticlesSlab_range(
                    rank, i, k_start_global, k_extent,
                    slab1, slab2, slab3, slab4, *param
                );
                processed++;
                printf("  Rank %d, i=%d: OK (X=[%d,%d))\n",
                       rank, i, k_start_global, k_start_global + k_extent);
            } catch (const std::exception& e) {
                fprintf(stderr, "ERROR: WriteParticlesSlab_range failed for rank %d, i=%d: %s\n",
                        rank, i, e.what());
            } catch (...) {
                fprintf(stderr, "ERROR: Unknown exception in WriteParticlesSlab_range for rank %d, i=%d\n",
                        rank, i);
            }
        }
    }
    
    printf("\n================================================================================\n");
    printf("COMPLETE: Processed %d rank-i pairs\n", processed);
    printf("Output files written to: %s\n", param->output_dir.c_str());
    printf("  Files: ic_rank*_i*_x*_* (per-rank particle data)\n");
    if (param->qdensity) {
        printf("  Density files: dens_rank*_i*_x*_* (per-rank density data)\n");
    }
    printf("================================================================================\n");
    
    // Cleanup
    TeardownOutput();
    delete param;
    
    return 0;
}

