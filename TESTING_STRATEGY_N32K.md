# Testing Strategy: Options A vs B for N=2048 (Medium-Scale Testing)

## Overview

This document provides a comprehensive testing strategy for comparing Options A and B with N=2048. The results will help you decide which option to use for production runs with N=32,000 (where Option B benefits are even larger).

## Recommended Testing Approach

### Phase 1: Small-Scale Validation (N=256, 512)

**Purpose**: Verify both options work correctly before scaling up.

1. **Quick Correctness Test**:
   ```bash
   # Test Option A
   make clean && make CFLAGS="-DUSE_PARTICLE_OUTPUT_OPTION_B=0"
   mpirun -np 16 ./main 256 examples/param_N256.par
   
   # Test Option B
   make clean && make CFLAGS="-DUSE_PARTICLE_OUTPUT_OPTION_B=1"
   mpirun -np 16 ./main 256 examples/param_N256.par
   
   # Compare outputs
   ./compare_option_outputs.sh output output --verify-data
   ```

2. **Verify outputs match**: Both should produce identical particle files.

### Phase 2: Medium-Scale Performance (N=2048, 4096)

**Purpose**: Measure performance differences before full-scale test.

1. **Run both options** with timing:
   ```bash
   ./test_options_A_vs_B_N32K.sh examples/param_N2048.par 256 16
   ```

2. **Key metrics to monitor**:
   - Total execution time
   - Memory usage per rank
   - Transpose time (Option A only)
   - I/O time

### Phase 3: Production-Scale Test (N=32,000) - Optional

**Purpose**: Final validation for production deployment (if needed). For N=32K, Option B benefits are even larger than N=2048.

## Testing Scripts

### 1. Automated Comparison Script

**File**: `test_options_A_vs_B_N32K.sh`

**Usage**:
```bash
./test_options_A_vs_B_N32K.sh [param_file] [num_ranks] [num_nodes]
```

**What it does**:
- Compiles both Option A and Option B
- Runs each with timing and memory monitoring
- Compares outputs
- Generates summary report

**Example**:
```bash
# For N=2048 with 16×16 ranks on 16 nodes
./test_options_A_vs_B_N32K.sh examples/param_N2048.par 256 16
```

### 2. Output Comparison Script

**File**: `compare_option_outputs.sh`

**Usage**:
```bash
./compare_option_outputs.sh output_option_a output_option_b [--verify-data]
```

**What it does**:
- Compares file counts
- Compares file sizes
- Optionally verifies data matches (binary comparison)

## SLURM Job Submission (Recommended for N=32K)

For N=32,000, you'll likely need to submit SLURM jobs. Here's a template:

### Option A Test Job (N=2048)

```bash
#!/bin/bash
#SBATCH --job-name=test_option_a_N2048
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB              # Adjust based on Option A memory needs
#SBATCH --time=01:00:00
#SBATCH --output=option_a_%j.out
#SBATCH --error=option_a_%j.err

# Compile Option A
cd $SLURM_SUBMIT_DIR
make clean
make CFLAGS="-DUSE_PARTICLE_OUTPUT_OPTION_B=0"

# Run
srun ./main 2048 examples/param_N2048.par
```

### Option B Test Job (N=2048)

```bash
#!/bin/bash
#SBATCH --job-name=test_option_b_N2048
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB              # Less memory needed for Option B
#SBATCH --time=00:45:00
#SBATCH --output=option_b_%j.out
#SBATCH --error=option_b_%j.err

# Compile Option B
cd $SLURM_SUBMIT_DIR
make clean
make CFLAGS="-DUSE_PARTICLE_OUTPUT_OPTION_B=1"

# Run
srun ./main 2048 examples/param_N2048.par
```

## Key Metrics to Monitor

### 1. Memory Usage

**Option A (N=2048, 16×16 ranks)**:
- `local_z_slab`: ~0.05 GB per rank
- Transposed slabs: ~0.05 GB per rank (transient)
- **Peak**: ~0.1 GB per rank

**Option B (N=2048, 16×16 ranks)**:
- `local_z_slab`: ~0.05 GB per rank
- **Peak**: ~0.05 GB per rank

**For N=32K (production scale)**:
- Option A: ~1.7 GB per rank
- Option B: ~0.85 GB per rank

**Measurement**:
```bash
# In SLURM job, check memory usage
sstat -j $SLURM_JOB_ID --format=MaxRSS,MaxVMSize
```

### 2. Execution Time

**Components to measure**:
- Total time
- Transpose time (Option A only) - look for "transposed and called" messages
- I/O time (particle writing)

**Expected for N=2048**:
- Option A: ~20s per rank (includes ~8s transpose)
- Option B: ~12s per rank (no transpose)

**For N=32K (production scale)**:
- Option A: ~135s per rank (includes ~60s transpose)
- Option B: ~75s per rank (no transpose)

### 3. Correctness

**Verification steps**:
1. **File counts**: Both should produce same number of files
2. **File sizes**: Should match exactly
3. **Particle data**: Sample files should be identical (binary comparison)
4. **Global k indices**: Verify `k_global = k_start_global + k_local` is correct

## Recommended Test Sequence

### Step 1: Quick Validation (N=256, 1 hour)

```bash
# Test both options
./test_options_A_vs_B_N32K.sh examples/param_N256.par 16 1

# Verify outputs match
./compare_option_outputs.sh output_option_a output_option_b --verify-data
```

**Success criteria**:
- Both complete without errors
- Output files match (count, size, data)

### Step 2: Medium-Scale Performance (N=2048, 4 hours)

```bash
# Submit SLURM jobs for both options
sbatch test_option_a_N2048.slurm
sbatch test_option_b_N2048.slurm

# After completion, compare
./compare_option_outputs.sh output_option_a output_option_b
```

**Success criteria**:
- Option B is faster (no transpose overhead)
- Option B uses less memory
- Outputs are identical

### Step 3: Production-Scale Test (N=32K, 1-2 days) - Optional

**Note**: If N=2048 tests confirm Option B benefits, you can proceed to production with N=32K. Option B benefits scale with N, so if it's better for N=2048, it will be even better for N=32K.

```bash
# Submit production-scale jobs (if needed for final validation)
sbatch test_option_a_N32K.slurm  # Edit to use N=32000
sbatch test_option_b_N32K.slurm  # Edit to use N=32000

# After completion, analyze results
./analyze_production_results.sh
```

**Success criteria**:
- Option B saves significant time (~60s per Z-slab)
- Option B saves significant memory (~0.85 GB per rank)
- Outputs are identical
- No memory errors or OOM kills

## Analysis Script

Create `analyze_production_results.sh` to:

1. **Extract timing from logs**:
   ```bash
   grep "Total time" option_a_*.out
   grep "Total time" option_b_*.out
   ```

2. **Extract memory from SLURM**:
   ```bash
   sacct -j $JOB_ID_A --format=MaxRSS,Elapsed
   sacct -j $JOB_ID_B --format=MaxRSS,Elapsed
   ```

3. **Calculate savings**:
   - Time saved = Option A time - Option B time
   - Memory saved = Option A memory - Option B memory
   - CPU-hours saved = (time saved × num_ranks) / 3600

## Decision Matrix

### Choose Option B if:
- ✅ N ≥ 8192 (your case: N=32K)
- ✅ Memory per node is constrained
- ✅ Performance is critical
- ✅ Test results show significant savings

### Choose Option A if:
- ✅ N < 8192
- ✅ Memory is abundant
- ✅ Code simplicity is priority
- ✅ Test results show minimal difference

## Expected Results for N=2048

Based on calculations:

| Metric | Option A | Option B | Improvement |
|--------|----------|----------|-------------|
| **Memory/rank** | ~0.1 GB | ~0.05 GB | **2× less** |
| **Transpose time** | ~8s/Z-slab | 0s | **8s saved** |
| **Total time/rank** | ~20s | ~12s | **40% faster** |
| **CPU-hours (256 ranks)** | ~1.4 h | ~0.85 h | **~0.6 h saved** |

## Expected Results for N=32K (Production Scale)

| Metric | Option A | Option B | Improvement |
|--------|----------|----------|-------------|
| **Memory/rank** | ~1.7 GB | ~0.85 GB | **2× less** |
| **Transpose time** | ~60s/Z-slab | 0s | **60s saved** |
| **Total time/rank** | ~135s | ~75s | **45% faster** |
| **CPU-hours (6561 ranks)** | 246 h | 137 h | **109 h saved** |

## Troubleshooting

### Issue: "Out of memory" errors with Option A
**Solution**: Switch to Option B (saves 0.85 GB per rank)

### Issue: Outputs don't match
**Check**:
1. Same parameter file used?
2. Same random seed?
3. Same number of ranks?
4. Grid decomposition matches?

### Issue: Option B slower than expected
**Check**:
1. Cache locality issues?
2. I/O bottleneck?
3. Network contention?

## Next Steps After Testing

1. **If Option B wins** (expected for N=32K):
   - Update `config.h`: `#define USE_PARTICLE_OUTPUT_OPTION_B 1`
   - Recompile for production
   - Document the decision

2. **If Option A wins** (unlikely for N=32K):
   - Keep default setting
   - Document why (e.g., memory not an issue)

3. **Production deployment**:
   - Use tested option
   - Monitor first production run
   - Keep Option B code available for future scaling

## Summary

For N=2048, **Option B is recommended** based on:
- **Memory savings**: ~0.05 GB per rank × 256 ranks = ~13 GB total saved
- **Time savings**: ~8s per Z-slab × 128 Z-slabs = ~17 minutes saved per run
- **Scalability**: Benefits increase with larger N

For N=32,000 (production), **Option B is strongly recommended**:
- **Memory savings**: 0.85 GB per rank × 6561 ranks = 5.6 TB total saved
- **Time savings**: ~60s per Z-slab × 81 Z-slabs = ~81 minutes saved per run
- **Scalability**: Better for production-scale runs

**Testing with N=2048 validates these benefits and ensures correctness. The benefits scale with N, so if Option B is better for N=2048, it will be even better for N=32K.**

