# PLT-Disabled Test Directory

Tests with PLT (Particle Linear Theory) disabled to isolate the Z-correlation issue.

## Files

- `N16_assembly_no_plt.pbs`: PBS script for running N=16 test with PLT disabled
- `param_N16_no_plt.par`: Parameter file with `ZD_qPLT = 0`
- `README.md`: This file

## Purpose

Test if the Z-correlation issue persists when PLT is disabled. This will help determine if:
1. The issue is specific to PLT mode
2. The issue is in the non-PLT displacement computation
3. There's a coordinate system convention difference

## Usage

```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly_no_plt
qsub N16_assembly_no_plt.pbs
```

## Debug Output

The compilation includes `-DDEBUG_EIGENVECTOR=1` which will print eigenvector components for test coordinates (x,y,z <= 2) when PLT is enabled. Since PLT is disabled in this test, the debug output will show the k-vector components instead.

## Running zeldovich-PLT with PLT Disabled

For a fair comparison, you should also run zeldovich-PLT with PLT disabled:

```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly_no_plt
./run_zeldovich_no_plt.sh
```

This will:
- Run zeldovich-PLT with the same parameter file (PLT disabled)
- Output to `/home/helenshao/InitialConditions/zeldovich-PLT/output_no_plt/`
- Generate reference output for comparison

## Comparison

After running both tests, compare the outputs:

1. **PLT Disabled Comparison:**
   ```bash
   # Compare hermitian (PLT disabled) vs zeldovich (PLT disabled)
   python3 ../test_assembly/compare_ic_files.py \
       ../test_assembly_no_plt/particle_ics \
       /home/helenshao/InitialConditions/zeldovich-PLT/output_no_plt
   ```

2. **PLT Enabled Comparison:**
   ```bash
   # Compare hermitian (PLT enabled) vs zeldovich (PLT enabled)
   python3 ../test_assembly/compare_ic_files.py \
       ../test_assembly/particle_ics \
       /home/helenshao/InitialConditions/zeldovich-PLT/output
   ```

This will help determine if the Z-correlation issue:
- Exists in both PLT and non-PLT modes
- Is specific to PLT mode
- Is specific to non-PLT mode

