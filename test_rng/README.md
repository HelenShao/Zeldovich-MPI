# RNG Comparison Test

Verify that both codes generate the same random Gaussian numbers (D values) before any FFT operations. This helps ensure RNG consistency and identify any differences in random number generation.

## Test Configuration

- **Grid size**: N=16
- **Parallelization**: Serial (1 MPI rank, 1 OMP thread)
- **PLT**: Disabled (ZD_qPLT = 0)
- **Random seed**: 4 (matches in both codes)
- **Debug output**: [RNG-DEBUG] lines for coordinates x,y,z <= 10

## Files

- `N16_rng_test.pbs`: PBS script to run hermitian_3d_matrix (1 rank, 1 thread)
- `param_N16_rng.par`: Parameter file with PLT disabled
- `run_zeldovich_rng.sh`: Script to run zeldovich-PLT (serial)
- `workflow.txt`: Step-by-step instructions for running the test
- `rng_logs/`: Directory containing filtered RNG debug output
  - `hermitian_rng_debug_filtered.txt`: RNG debug from hermitian code
  - `zeldovich_rng_debug_filtered.txt`: RNG debug from zeldovich code

## Expected Output

Both codes should print [RNG-DEBUG] lines in the format:
```
[RNG-DEBUG] N=16 Y=0 (x,z)=(0,0): k=(0,0,0) k2=0.000000 | D=(0.000000e+00,0.000000e+00) F=(...) G=(...) H=(...)
```

For N=16, this will print approximately:
- 1,331 coordinates per Y value (x,y,z <= 10)
- 16 Y values total
- ~21,296 total lines (some may be skipped for DC/Nyquist modes)

## Running the Test

See `workflow.txt` for detailed instructions.

Quick start:
1. `qsub N16_rng_test.pbs` (run hermitian code)
2. `./run_zeldovich_rng.sh` (run zeldovich code)
3. Compare: `diff rng_logs/hermitian_rng_debug_filtered.txt rng_logs/zeldovich_rng_debug_filtered.txt`

## What to Compare

1. **D values**: Should match exactly (same random Gaussian draws)
2. **F, G, H values**: Should match if D matches and computation is identical
3. **k, k2 values**: Should match (same k-vector computation)
4. **Line count**: Should be approximately the same (may differ slightly due to DC/Nyquist handling)

## Notes

- Both codes use the same random seed (4) and same parameter file
- PLT is disabled to simplify comparison
- Serial execution (1 rank, 1 thread) eliminates MPI/OpenMP differences
- RNG debug output is printed before any FFT operations

