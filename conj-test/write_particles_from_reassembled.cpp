// ====================================================================================
// WRITE PARTICLES FROM REASSEMBLED Z-SLABS
// ====================================================================================
// Reads reassembled Z-slabs and calls WriteParticlesSlab to write particle data
// This program:
// 1. Reads command-line arguments
// 2. Loads or creates simulation parameters
// 3. Creates a `BlockArray` object (required by `WriteParticlesSlab`)
// 4. Initializes output buffers
// 5. For each Z-slab:
//    - Reads the reassembled binary file
//    - Extracts pointers to the 4 arrays
//    - Calls `WriteParticlesSlab` to write particle data
//
// Example:
//   ./write_particles_from_reassembled reassembled_z_slabs 16 param_N16.par 0 16
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
#include <block_array.h>
#include <parameters.h>
#include <zeldovich.h>
#include <omp.h>

namespace fs = std::filesystem;
using Complx = std::complex<double>;

// ====================================================================================
// READ REASSEMBLED Z-SLAB
// ====================================================================================

std::vector<Complx> read_reassembled_z_slab(const std::string& filename, int N, int narray) {
    std::vector<Complx> full_slab;
    
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open file: %s\n", filename.c_str());
        return full_slab;
    }
    
    // Get file size
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    
    // Expected size: narray * N * N * 16 bytes (double precision)
    size_t expected_size = narray * N * N * 16;
    if (file_size != (long)expected_size) {
        fprintf(stderr, "ERROR: File size mismatch for %s: expected %zu, got %ld\nCheck precision!",
                filename.c_str(), expected_size, file_size);
        fclose(fp);
        return full_slab;
    }
    
    // Allocate and read
    full_slab.resize(narray * N * N);
    size_t read_size = fread(full_slab.data(), 16, narray * N * N, fp);
    fclose(fp);
    
    if (read_size != narray * N * N) {
        fprintf(stderr, "ERROR: Read %zu elements, expected %d\n", 
                read_size, narray * N * N);
        full_slab.clear();
    }
    
    return full_slab;
}

// ====================================================================================
// MAIN PROGRAM
// ====================================================================================

int main(int argc, char* argv[]) {
    // Force single OpenMP thread to avoid segfault in zeldovich-PLT initialization
    omp_set_num_threads(1);
    
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <reassembled_dir> <N> [param_file] [z_start] [z_end]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  reassembled_dir: Directory containing z*_slab_N*_full.bin files\n");
        fprintf(stderr, "  N: Resolution\n");
        fprintf(stderr, "  param_file: Optional parameter file (default: create minimal)\n");
        fprintf(stderr, "  z_start: First Z-slab to process (default: 0)\n");
        fprintf(stderr, "  z_end: Last Z-slab to process (default: N)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Example:\n");
        fprintf(stderr, "  %s reassembled_z_slabs 16 param_N16.par 0 16\n", argv[0]);
        return 1;
    }
    
    std::string reassembled_dir = argv[1];
    int N = std::stoi(argv[2]);
    int narray = 4;
    
    std::string param_file = "";
    int z_start = 0;
    int z_end = N;
    
    if (argc >= 4) {
        param_file = argv[3];
    }
    if (argc >= 5) {
        z_start = std::stoi(argv[4]);
    }
    if (argc >= 6) {
        z_end = std::stoi(argv[5]);
    }
    
    printf("================================================================================\n");
    printf("WRITE PARTICLES FROM REASSEMBLED Z-SLABS\n");
    printf("================================================================================\n");
    printf("Reassembled directory: %s\n", reassembled_dir.c_str());
    printf("Grid size N: %d\n", N);
    printf("Number of arrays: %d\n", narray);
    printf("Z-slab range: [%d, %d)\n", z_start, z_end);
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
        std::string tmp_param_file = "/tmp/write_particles_params.par";
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
            fprintf(tmp_fp, "InitialConditionsDirectory = \"./output\"\n");
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
    
    // Create BlockArray
    printf("Creating BlockArray...\n");
    if (N % 2 != 0) {
        fprintf(stderr, "ERROR: N=%d must be even for BlockArray (numblock=2)\n", N);
        delete param;
        return 1;
    }
    
    BlockArray* array = nullptr;
    try {
        fs::path tmp_dir = "/tmp/write_particles_blockarray";
        fs::create_directories(tmp_dir);
        array = new BlockArray(N, 2, narray, tmp_dir, 0, 1, -1);
        printf("  Created BlockArray (ppd=%ld, numblock=%d)\n", array->ppd, array->numblock);
    } catch (const std::exception& e) {
        fprintf(stderr, "ERROR: Failed to create BlockArray: %s\n", e.what());
        delete param;
        return 1;
    } catch (...) {
        fprintf(stderr, "ERROR: Unknown exception creating BlockArray\n");
        delete param;
        return 1;
    }
    
    // Initialize output buffers (required by WriteParticlesSlab)
    printf("Initializing output...\n");
    SetupOutputDir(*param);
    double buffer_size = InitOutputBuffers(*param);
    printf("  Output directory: %s\n", param->output_dir.c_str());
    printf("  Buffer size: %.3f GiB\n", buffer_size);
    
    // WriteParticlesSlab writes to internal buffers, FILE* can be NULL
    FILE* output_fp = NULL;
    
    printf("\nProcessing Z-slabs...\n");
    
    // Process each Z-slab
    int processed = 0;
    for (int z = z_start; z < z_end; z++) {
        // Read reassembled Z-slab
        std::string filename = reassembled_dir + "/z" + std::to_string(z) + "_slab_N" + std::to_string(N) + "_full.bin";
        std::vector<Complx> full_slab = read_reassembled_z_slab(filename, N, narray);
        
        if (full_slab.empty()) {
            printf("  Z=%d: Skipping (file not found or error)\n", z);
            continue;
        }
        
        // Verify we have the right amount of data
        if (full_slab.size() != (size_t)(narray * N * N)) {
            fprintf(stderr, "ERROR: Z=%d wrong size: %zu != %d\n", z, full_slab.size(), narray * N * N);
            continue;
        }
        
        // Create 2D slab pointers for WriteParticlesSlab
        // WriteParticlesSlab expects: slab[x + ppd * y] (row-major, Y varies slowly)
        Complx* slab1 = &full_slab[0 * N * N];  // Array 0: D + i*F
        Complx* slab2 = &full_slab[1 * N * N];  // Array 1: G + i*H
        Complx* slab3 = &full_slab[2 * N * N];  // Array 2: X-velocity
        Complx* slab4 = &full_slab[3 * N * N];  // Array 3: Y-velocity + Z-velocity
        
        // Call WriteParticlesSlab
        try {
            WriteParticlesSlab(output_fp, z, slab1, slab2, slab3, slab4, *array, *param);
            processed++;
            printf("  Z=%d: OK\n", z);
        } catch (const std::exception& e) {
            fprintf(stderr, "ERROR: WriteParticlesSlab failed for Z=%d: %s\n", z, e.what());
        } catch (...) {
            fprintf(stderr, "ERROR: Unknown exception in WriteParticlesSlab for Z=%d\n", z);
        }
    }
    
    // Note: output_fp is NULL since WriteParticlesSlab uses internal buffers
    // Don't call fclose(NULL) - it causes undefined behavior
    
    printf("\n================================================================================\n");
    printf("COMPLETE: Processed %d Z-slabs\n", processed);
    printf("Output files written to: %s\n", param->output_dir.c_str());
    printf("  Files: ic_* (particle data)\n");
    if (param->qdensity) {
        printf("  Density file: %s\n", param->density_filename.c_str());
    }
    printf("================================================================================\n");
    
    // Cleanup
    TeardownOutput();
    delete array;
    delete param;
    
    return 0;
}
