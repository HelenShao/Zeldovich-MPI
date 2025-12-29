# Issue Summary: Poor Z-Correlation in Particle IC Output

## Problem Description

The `hermitian_3d_matrix_production` code generates particle initial conditions (ICs) that show **poor Z-correlation** (0.06) with the reference `zeldovich-PLT` code, while X and Y correlations are reasonable (0.81 and 0.96 respectively). This occurs when using **Path 2** particle IC writing workflow (generating `.bin` files, then reassembling them).

### Current Status

- **Large values issue**: RESOLVED (fixed by reverting eigenvector normalization)
- **Displacement ranges**: Now reasonable for all components (X, Y, Z)
- **Correlations with zeldovich-PLT**:
  - X correlation: ~0.81 (moderate)
  - Y correlation: ~0.96 (good)
  - Z correlation: ~0.06 (very poor) ← **REMAINING ISSUE**
- **File sizes**: `.bin` files are 4096 bytes (correct for single precision: 4 arrays × 16 Y × 8 X × 8 bytes)

## Test Configuration

### Current Setup
- **Grid size**: N=16
- **Precision**: Single precision (float)
- **PLT**: Enabled (`ZD_qPLT=1` in `param_N16.par`)
- **PLT rescaling**: Enabled (`ZD_qPLTrescale=1`)
- **Parameter file**: `examples/param_N16.par`
- **Random seed**: 4
- **Output mode**: Mode 1 (write `.bin` files for reassembly)

### Test Locations

1. **Main test directory**: `/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly/`
   - PBS script: `N16_assembly.pbs`
   - Output: `particle_ics/` (IC files), `bin_files/` (intermediate `.bin` files)
   - Logs: `N16_bin_generation.log`, `N16_reassembly.log`

2. **PLT-disabled test directory**: `/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly_no_plt/`
   - PBS script: `N16_assembly_no_plt.pbs`
   - Parameter file: `param_N16_no_plt.par` (`ZD_qPLT=0`)
   - Purpose: Isolate PLT vs non-PLT issues

3. **Comparison scripts**:
   - `compare_ic_files.py`: Compares particle values between hermitian and zeldovich-PLT outputs
   - `visualize_particle_comparison.py`: Generates correlation plots and histograms

4. **Reference output**: `/home/helenshao/InitialConditions/zeldovich-PLT/output/` (PLT enabled)
   - `/home/helenshao/InitialConditions/zeldovich-PLT/output_no_plt/` (PLT disabled)

## Tests Performed

### 1. Initial Comparison (PLT Enabled)
- **Result**: Large values in Z and Y, poor Z-correlation (0.06)
- **Finding**: Normalization mismatch identified (~250× amplitude difference)
- **Fix applied**: Corrected `factor` computation in `hermitian_generation.c` (changed from `rescale * fundamental * ik2` to `rescale / (k2 * fundamental)`)
- **Outcome**: Y and X correlations improved (0.96 and 0.81), but Z-correlation remained poor

### 2. Eigenvector Investigation
- **Hypothesis**: Z-correlation issue due to eigenvector component mapping
- **Test**: Extracted eigenvector components from F, G, H values in logs
- **Result**: Eigenvectors match expected values (e.g., k=(0,3,0): e.vec[0]=0, e.vec[1]=1.0, e.vec[2]=0)
- **Conclusion**: Eigenvector retrieval appears correct

### 3. PLT-Disabled Test
- **Purpose**: Isolate PLT-specific issues
- **Setup**: Created `test_assembly_no_plt/` with PLT disabled
- **Result**: Z-correlation issue persists in both PLT and non-PLT modes
- **Finding**: Issue is in fundamental displacement computation, not PLT-specific

### 4. Eigenvector Normalization Attempt
- **Hypothesis**: Missing normalization `norm = k2 / (k·ehat)` as in zeldovich-PLT
- **Change**: Added normalization in 3 locations in `hermitian_generation.c`
- **Result**: Introduced very large values (~2.4×10¹⁰) in Z and Y displacements
- **Action**: Reverted changes (code restored to original state)
- **Outcome**: Large values disappeared after revert, confirming eigenvector normalization was the cause

### 5. Precision Mismatch Investigation
- **Finding**: `.bin` files are single precision (4096 bytes correct)
- **Issue**: Reassembly tool might read with wrong precision
- **Fix**: Updated PBS script to compile reassembly tool with matching CFLAGS
- **Result**: Confirmed precision is consistent (not the source of issues)

## Current State

### Code Status
- **Eigenvector normalization**: Reverted (using `ehat` directly from `plt_get_eigenmode`)
- **Normalization factor**: Fixed (`rescale / (k2 * fundamental)`)
- **Precision**: Single precision (both main and reassembly tool)
- **PBS script**: Updated to ensure matching precision flags

### Known Issues
1. **Z-correlation**: Very poor (0.06) even after normalization factor fix - **PRIMARY REMAINING ISSUE**
2. **Array mapping**: Unclear if arrays are written in correct order (Array 0: D+i*F, Array 1: G+i*H, Array 2: X-vel, Array 3: Y+Z-vel)
3. **Eigenvector usage**: Using `ehat` directly from `plt_get_eigenmode` (without `k2/(k·ehat)` normalization) - this is correct for hermitian code, different from zeldovich-PLT

### Files Modified
- `src/generation/hermitian_generation.c`: Normalization factor fix (reverted eigenvector normalization)
- `test_assembly/N16_assembly.pbs`: Updated to ensure matching precision flags
- `src/write_particles_from_reassembled_mpi.cpp`: Precision-aware reading (BinComplx)

## Future Tests

### 1. Serial Version Test (Single MPI Rank, Single Thread)
**Purpose**: Eliminate MPI communication errors as a source of the problem

**Steps**:
```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production
# Create test directory
mkdir -p test_assembly_serial
cd test_assembly_serial

# Run with 1 rank, 1 thread
export OMP_NUM_THREADS=1
./hermitian_3d_matrix 16 ../examples/param_N16.par

# Compare with zeldovich-PLT serial output
# (zeldovich-PLT should also be run serially for fair comparison)
```

**Expected outcome**: If Z-correlation improves significantly, the issue may be MPI-related. If it remains poor, the issue is in the core computation logic (displacement computation or array mapping).

### 2. Small Grid Test with PLT Disabled (N=8)
**Purpose**: Verify random number generation and density field consistency

**Steps**:
```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production
mkdir -p test_assembly_N8_no_plt
cd test_assembly_N8_no_plt

# Create param file with N=8, PLT disabled
# Copy param_N16_no_plt.par and modify:
#   CPD = 8
#   NP = 512  # 8^3

# Run hermitian_3d_matrix
./hermitian_3d_matrix 8 param_N8_no_plt.par

# Add debug output to print density field:
#   - Before FFT (in hermitian_generation.c)
#   - After FFT (in main.cpp before writing)
#   - Compare with zeldovich-PLT density field
```

**Debug additions needed**:
- In `hermitian_generation.c`: Print density field D for first few k-vectors
- In `main.cpp`: Print density field after 1D FFT for first Z-slab
- Compare with zeldovich-PLT density output

**Expected outcome**: 
- If density fields match: Random number generation is correct, issue is in displacement computation
- If density fields differ: Random number generation or power spectrum normalization is wrong

### 3. Array Order Verification
**Purpose**: Verify arrays are written in correct order to `.bin` files

**Steps**:
- Add debug output in `main.cpp` when writing `.bin` files (Mode 1)
- Print first few values of each array (0, 1, 2, 3) before writing
- Compare with what `write_particles_from_reassembled_mpi.cpp` expects:
  - Array 0: D + i*F (density + i*X displacement)
  - Array 1: G + i*H (Y displacement + i*Z displacement)
  - Array 2: X velocity
  - Array 3: Y velocity + Z velocity

### 4. Direct Comparison with zeldovich-PLT Array Values
**Purpose**: Compare intermediate array values (F, G, H) before writing

**Steps**:
- Add debug output in `hermitian_generation.c` to print F, G, H values for specific k-vectors
- Add corresponding debug output in `zeldovich-PLT/src/zeldovich.cpp`
- Compare values directly (not just final particle positions)

### 5. Check Array Initialization
**Purpose**: Verify all arrays are properly initialized before use

**Steps**:
- Check if `local_z_slab` is properly zeroed before unpacking
- Verify arrays 2 and 3 (velocities) are computed correctly
- Check if uninitialized arrays could cause garbage values

## Key Files and Locations

### Source Code
- `src/generation/hermitian_generation.c`: Core displacement computation (F, G, H)
- `src/main.cpp`: Main loop, `.bin` file writing (Mode 1)
- `src/write_particles_from_reassembled_mpi.cpp`: Reassembly tool
- `src/output/output_new.cpp`: Particle writing logic (`WriteParticlesSlab_new`)

### Configuration
- `examples/param_N16.par`: PLT-enabled parameter file
- `test_assembly_no_plt/param_N16_no_plt.par`: PLT-disabled parameter file
- `src/config.h`: Compile-time flags
- `src/precision.h`: Precision selection (single/double)

### Test Scripts
- `test_assembly/N16_assembly.pbs`: Main PBS script (PLT enabled)
- `test_assembly_no_plt/N16_assembly_no_plt.pbs`: PLT-disabled PBS script
- `test_assembly/compare_ic_files.py`: Particle comparison script
- `test_assembly/visualize_particle_comparison.py`: Visualization script

### Output Directories
- `test_assembly/particle_ics/`: Final particle IC files
- `test_assembly/bin_files/`: Intermediate `.bin` files
- `test_assembly_no_plt/particle_ics/`: PLT-disabled particle ICs

## Debugging Tips

1. **Check `.bin` file contents directly**:
   ```python
   import struct
   with open('test_assembly/bin_files/rank_0/i0_slab_N16.bin', 'rb') as f:
       data = f.read(32)  # First 4 complex floats
       for i in range(4):
           re, im = struct.unpack('ff', data[i*8:(i+1)*8])
           print(f"Array 0, Y=0, X={i}: ({re:.6f}, {im:.6f})")
   ```

2. **Compare with zeldovich-PLT intermediate values**:
   - Add debug prints in both codes for same k-vectors
   - Compare F, G, H values before FFT

3. **Verify array mapping**:
   - Check `output_new.cpp` line 184-186: `pos[0] = imag(slab2)`, `pos[1] = real(slab2)`, `pos[2] = imag(slab1)`
   - Verify this matches expected mapping: Z=imag(slab2), Y=real(slab2), X=imag(slab1)

## Notes

- **Large values issue**: RESOLVED by reverting eigenvector normalization. The `k2/(k·ehat)` normalization used in zeldovich-PLT is NOT appropriate for hermitian_3d_matrix_production.
- **Normalization factor fix**: The `rescale / (k2 * fundamental)` fix improved X and Y correlations (0.81 and 0.96) but Z-correlation remains poor (0.06).
- **Z-correlation issue**: The poor Z-correlation persists in both PLT-enabled and PLT-disabled modes, suggesting it's not PLT-specific but rather in the fundamental displacement computation or array mapping.
- **Eigenvector usage**: The code correctly uses `ehat` directly from `plt_get_eigenmode` without additional normalization. The difference from zeldovich-PLT is intentional and correct.

