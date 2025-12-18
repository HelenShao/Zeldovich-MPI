// ====================================================================================
// REASSEMBLE Z-SLABS AND CALL WriteParticlesSlab
// ====================================================================================
// This program reads distributed rank files (rank_*/z*_slab_N*.bin) and reassembles
// them into full Z-slabs, then calls WriteParticlesSlab to write particle data.
//
// Usage:
//   ./reassemble_z_slabs <output_dir> <N> [param_file] [z_start] [z_end]
//
// Example:
//   ./reassemble_z_slabs ./output 256 param.par 0 256
// ====================================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <vector>
#include <string>
#include <glob.h>
#include <sys/stat.h>
#include <filesystem>
#include <iostream>

// Include zeldovich-PLT headers
#include <output.h>
#include <block_array.h>
#include <parameters.h>
#include <zeldovich.h>

namespace fs = std::filesystem;
using Complx = std::complex<double>;

// ====================================================================================
// FILE READING UTILITIES
// ====================================================================================

struct RankFile {
    int rank;
    int x_count;
    std::vector<Complx> data;  // [narray][N][x_count] in memory
};

// Read a single rank file
// Format: [Array][Y][X] as fftw_complex_t (8 or 16 bytes depending on precision)
RankFile read_rank_file(const std::string& filename, int N, int narray) {
    RankFile result;
    result.rank = -1;
    
    // Parse rank from filename: rank_*/z*_slab_N*.bin
    size_t rank_pos = filename.find("rank_");
    if (rank_pos == std::string::npos) {
        fprintf(stderr, "ERROR: Cannot parse rank from filename: %s\n", filename.c_str());
        return result;
    }
    
    size_t slash_pos = filename.find('/', rank_pos);
    if (slash_pos == std::string::npos) {
        fprintf(stderr, "ERROR: Cannot parse rank from filename: %s\n", filename.c_str());
        return result;
    }
    
    std::string rank_str = filename.substr(rank_pos + 5, slash_pos - rank_pos - 5);
    result.rank = std::stoi(rank_str);
    
    // Open file
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open file: %s\n", filename.c_str());
        return result;
    }
    
    // Get file size
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    
    // Determine precision from file size
    // Expected size: narray * N * x_count * bytes_per_complex
    // Try both single (8 bytes) and double (16 bytes) precision
    int bytes_per_complex = 16;  // Default to double precision (more common)
    int x_count_single = file_size / (narray * N * 8);
    int x_count_double = file_size / (narray * N * 16);
    
    // Check which one gives exact match
    // Prefer double precision if both match (more common for zeldovich-PLT compatibility)
    bool single_match = (file_size == narray * N * x_count_single * 8);
    bool double_match = (file_size == narray * N * x_count_double * 16);
    
    if (double_match) {
        bytes_per_complex = 16;
        result.x_count = x_count_double;
    } else if (single_match) {
        bytes_per_complex = 8;
        result.x_count = x_count_single;
    } else {
        fprintf(stderr, "WARNING: File size %ld does not match expected format for %s\n", 
                file_size, filename.c_str());
        fprintf(stderr, "  Trying single precision: x_count=%d, remainder=%ld\n",
                x_count_single, file_size % (narray * N * 8));
        fprintf(stderr, "  Trying double precision: x_count=%d, remainder=%ld\n",
                x_count_double, file_size % (narray * N * 16));
        // Default to double precision
        bytes_per_complex = 16;
        result.x_count = x_count_double;
    }
    
    // Allocate buffer for reading
    size_t num_elements = narray * N * result.x_count;
    std::vector<char> buffer(file_size);
    
    // Read file
    size_t read_size = fread(buffer.data(), 1, file_size, fp);
    fclose(fp);
    
    if (read_size != (size_t)file_size) {
        fprintf(stderr, "ERROR: Read %zu bytes, expected %ld bytes from %s\n",
                read_size, file_size, filename.c_str());
        result.rank = -1;
        return result;
    }
    
    // Convert to Complx (std::complex<double>)
    result.data.resize(num_elements);
    
    if (bytes_per_complex == 8) {
        // Single precision: float[2]
        float* fdata = (float*)buffer.data();
        for (size_t i = 0; i < num_elements; i++) {
            result.data[i] = Complx((double)fdata[2*i], (double)fdata[2*i+1]);
        }
    } else {
        // Double precision: double[2]
        double* ddata = (double*)buffer.data();
        for (size_t i = 0; i < num_elements; i++) {
            result.data[i] = Complx(ddata[2*i], ddata[2*i+1]);
        }
    }
    
    return result;
}

// Reassemble full Z-slab from all rank files
// Returns: [narray][N][N] array, or empty if error
std::vector<Complx> reassemble_z_slab(int z, const std::string& output_dir, int N, int narray) {
    // Find all files for this Z-slab
    std::string pattern = output_dir + "/rank_*/z" + std::to_string(z) + "_slab_N" + std::to_string(N) + ".bin";
    
    glob_t glob_result;
    glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    
    if (glob_result.gl_pathc == 0) {
        fprintf(stderr, "WARNING: No files found for pattern: %s\n", pattern.c_str());
        globfree(&glob_result);
        return std::vector<Complx>();
    }
    
    // Read all rank files
    std::vector<RankFile> rank_files;
    for (size_t i = 0; i < glob_result.gl_pathc; i++) {
        RankFile rf = read_rank_file(glob_result.gl_pathv[i], N, narray);
        if (rf.rank >= 0) {
            rank_files.push_back(rf);
        }
    }
    globfree(&glob_result);
    
    if (rank_files.empty()) {
        fprintf(stderr, "ERROR: No valid rank files found for Z=%d\n", z);
        return std::vector<Complx>();
    }
    
    // Sort by rank number (assuming contiguous X-distribution)
    std::sort(rank_files.begin(), rank_files.end(),
              [](const RankFile& a, const RankFile& b) { return a.rank < b.rank; });
    
    // Verify total X count
    int total_x_count = 0;
    for (const auto& rf : rank_files) {
        total_x_count += rf.x_count;
    }
    
    if (total_x_count != N) {
        fprintf(stderr, "WARNING: Z=%d total x_count=%d != N=%d\n", z, total_x_count, N);
        if (abs(total_x_count - N) > 2) {
            fprintf(stderr, "ERROR: X count mismatch too large, aborting\n");
            return std::vector<Complx>();
        }
    }
    
    // Allocate full slab: [narray][N][N]
    std::vector<Complx> full_slab(narray * N * N);
    
    // Reassemble: copy data from each rank file
    int x_start = 0;
    for (const auto& rf : rank_files) {
        int x_end = std::min(x_start + rf.x_count, N);
        int actual_x_count = x_end - x_start;
        
        if (actual_x_count > 0 && x_start < N) {
            // Copy data: [narray][N][x_count] -> [narray][N][N]
            // File format: [Array][Y][X] (row-major)
            // Memory layout: data[array_idx * N * x_count + y * x_count + x_idx]
            for (int array_idx = 0; array_idx < narray; array_idx++) {
                for (int y = 0; y < N; y++) {
                    for (int x_idx = 0; x_idx < actual_x_count; x_idx++) {
                        int src_idx = array_idx * N * rf.x_count + y * rf.x_count + x_idx;
                        int dst_idx = array_idx * N * N + y * N + (x_start + x_idx);
                        full_slab[dst_idx] = rf.data[src_idx];
                    }
                }
            }
        }
        x_start = x_end;
    }
    
    return full_slab;
}

// ====================================================================================
// MINIMAL BLOCKARRAY WRAPPER
// ====================================================================================

// Minimal BlockArray wrapper for WriteParticlesSlab
// We only need ppd, other fields can be dummy values
class MinimalBlockArray {
public:
    int64_t ppd;
    
    MinimalBlockArray(int64_t _ppd) : ppd(_ppd) {}
};

// ====================================================================================
// MAIN PROGRAM
// ====================================================================================

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <output_dir> <N> [param_file] [z_start] [z_end]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  output_dir: Directory containing rank_*/z*_slab_N*.bin files\n");
        fprintf(stderr, "  N: Grid size (e.g., 256, 512)\n");
        fprintf(stderr, "  param_file: Optional parameter file for zeldovich-PLT (default: create minimal)\n");
        fprintf(stderr, "  z_start: First Z-slab to process (default: 0)\n");
        fprintf(stderr, "  z_end: Last Z-slab to process (default: N)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Example:\n");
        fprintf(stderr, "  %s ./output 256 param.par 0 256\n", argv[0]);
        return 1;
    }
    
    std::string output_dir = argv[1];
    int N = std::stoi(argv[2]);
    int narray = 4;  // Always 4 arrays
    
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
    printf("REASSEMBLE Z-SLABS AND CALL WriteParticlesSlab\n");
    printf("================================================================================\n");
    printf("Output directory: %s\n", output_dir.c_str());
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
        // We need to create a minimal parameter file or use defaults
        // For now, we'll create a temporary parameter file with minimal settings
        std::string tmp_param_file = "/tmp/reassemble_params.par";
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
            fprintf(tmp_fp, "ICFormat = \"RV\"\n");
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
    
    // Verify ppd matches N
    if (param->ppd != N) {
        fprintf(stderr, "WARNING: Parameter ppd=%ld != N=%d, using N\n", param->ppd, N);
        // We can't modify ppd directly, but WriteParticlesSlab uses array.ppd
    }
    
    // Create minimal BlockArray
    // Note: BlockArray constructor requires numblock to be even and ppd divisible by numblock
    // We use numblock=2 (minimum even value) and ensure N is even
    if (N % 2 != 0) {
        fprintf(stderr, "ERROR: N=%d must be even for BlockArray (numblock=2)\n", N);
        delete param;
        return 1;
    }
    
    BlockArray* array = nullptr;
    try {
        // BlockArray constructor: BlockArray(int _ppd, int _numblock, int _narray, 
        //                                    const fs::path &_dir, int _ramdisk, 
        //                                    int _quickdelete, int _part)
        // numblock=2 (minimum even value), part=-1 (both parts, but we won't use them)
        fs::path tmp_dir = "/tmp/reassemble_blockarray";
        fs::create_directories(tmp_dir);
        array = new BlockArray(N, 2, narray, tmp_dir, 0, 1, -1);
        printf("Created BlockArray (ppd=%ld, numblock=%d)\n", array->ppd, array->numblock);
    } catch (const std::exception& e) {
        fprintf(stderr, "ERROR: Failed to create BlockArray: %s\n", e.what());
        delete param;
        return 1;
    }
    
    // Open output file for WriteParticlesSlab
    // WriteParticlesSlab writes to FILE* output (can be NULL for some formats)
    // For now, we'll use stdout or create an output file
    std::string output_filename = output_dir + "/particles_output.bin";
    FILE* output_fp = fopen(output_filename.c_str(), "wb");
    if (!output_fp) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", output_filename.c_str());
        delete array;
        delete param;
        return 1;
    }
    
    printf("\nProcessing Z-slabs...\n");
    
    // Process each Z-slab
    int processed = 0;
    for (int z = z_start; z < z_end; z++) {
        // Reassemble full Z-slab
        std::vector<Complx> full_slab = reassemble_z_slab(z, output_dir, N, narray);
        
        if (full_slab.empty()) {
            if (z % 10 == 0) {  // Print every 10th missing slab
                printf("  Z=%d: Skipping (files not found or error)\n", z);
            }
            continue;
        }
        
        // Create 2D slab pointers for WriteParticlesSlab
        // Layout: [narray][N][N] -> each array is [N][N] contiguous
        // WriteParticlesSlab expects: slab[x + ppd * y] (row-major, Y varies slowly)
        // Our layout: full_slab[array_idx * N * N + y * N + x]
        // This matches! (x is stride-1, y is stride-N, same as WriteParticlesSlab expects)
        Complx* slab1 = &full_slab[0 * N * N];  // Array 0: D + i*F
        Complx* slab2 = &full_slab[1 * N * N];  // Array 1: G + i*H
        Complx* slab3 = &full_slab[2 * N * N];  // Array 2: X-velocity
        Complx* slab4 = &full_slab[3 * N * N];  // Array 3: Y-velocity + Z-velocity
        
        // Call WriteParticlesSlab
        try {
            WriteParticlesSlab(output_fp, z, slab1, slab2, slab3, slab4, *array, *param);
            processed++;
            
            if (z % 10 == 0 || z == z_end - 1) {
                printf("  Z=%d: Processed successfully\n", z);
            }
        } catch (const std::exception& e) {
            fprintf(stderr, "ERROR: WriteParticlesSlab failed for Z=%d: %s\n", z, e.what());
        }
    }
    
    fclose(output_fp);
    
    printf("\n================================================================================\n");
    printf("COMPLETE: Processed %d Z-slabs\n", processed);
    printf("Output file: %s\n", output_filename.c_str());
    printf("================================================================================\n");
    
    // Cleanup
    delete array;
    delete param;
    
    return 0;
}
