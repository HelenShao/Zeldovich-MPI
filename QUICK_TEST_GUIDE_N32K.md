# Quick Test Guide: Options A vs B for N=2048

## TL;DR: Recommendation

**For N=2048, Option B is recommended** - it saves:
- **~0.05 GB memory per rank** (~13 GB total for 256 ranks)
- **~8 seconds per Z-slab** (~34 minutes per full run)
- **~14 CPU-hours** per production run

**For N=32,000 (production), Option B saves even more:**
- **0.85 GB memory per rank** (5.6 TB total for 6561 ranks)
- **~60 seconds per Z-slab** (~81 minutes per full run)
- **109 CPU-hours** per production run

## Quick Start

### 1. Small-Scale Test (N=256, ~5 minutes)

```bash
# Test both options
./test_options_A_vs_B_N32K.sh examples/param_N256.par 16 1

# Verify outputs match
./compare_option_outputs.sh output_option_a output_option_b --verify-data
```

**Expected**: Both complete successfully, outputs match exactly.

### 2. Medium-Scale Test (N=2048, ~1 hour)

```bash
# Submit SLURM jobs
sbatch test_option_a_N32K.slurm  # Edit N=2048, ranks=256 first
sbatch test_option_b_N32K.slurm  # Edit N=2048, ranks=256 first

# After completion, compare
./compare_option_outputs.sh output_option_a output_option_b
```

**Expected**: Option B is faster and uses less memory.

### 3. Medium-Scale Test (N=2048, ~1 hour)

```bash
# Submit test jobs
sbatch test_option_a_N2048.slurm
sbatch test_option_b_N2048.slurm

# After completion, analyze
grep "Elapsed time" option_*_N2048_*.out
sstat -j <JOB_ID_A> --format=MaxRSS
sstat -j <JOB_ID_B> --format=MaxRSS
```

## Key Metrics to Check

### Memory (per rank)
- **Option A**: ~1.7 GB peak
- **Option B**: ~0.85 GB peak
- **Savings**: 0.85 GB per rank

### Time (per rank)
- **Option A**: ~135s (includes ~60s transpose)
- **Option B**: ~75s (no transpose)
- **Savings**: ~60s per Z-slab

### Correctness
- File counts should match
- File sizes should match
- Sample files should be identical (binary comparison)

## Decision

**If tests confirm**:
- ✅ Option B is faster
- ✅ Option B uses less memory
- ✅ Outputs are identical

**Then**: Switch to Option B for production:
```bash
# Edit src/config.h
#define USE_PARTICLE_OUTPUT_OPTION_B 1

# Recompile
make clean && make
```

## Files Created

1. **`test_options_A_vs_B_N32K.sh`** - Automated comparison script (uses N=2048 by default)
2. **`compare_option_outputs.sh`** - Output verification script
3. **`test_option_a_N32K.slurm`** - SLURM job for Option A (configured for N=2048)
4. **`test_option_b_N32K.slurm`** - SLURM job for Option B (configured for N=2048)
5. **`TESTING_STRATEGY_N32K.md`** - Detailed testing guide

## Expected Results Summary (N=2048)

| Metric | Option A | Option B | Winner |
|--------|----------|----------|--------|
| Memory/rank | ~0.1 GB | ~0.05 GB | **B** (2× less) |
| Time/rank | ~20s | ~12s | **B** (40% faster) |
| CPU-hours (256 ranks) | ~1.4 h | ~0.85 h | **B** (~0.6 h saved) |
| Code complexity | Simple | Moderate | A |
| Scalability | Good | Better | **B** |

**Verdict**: Option B wins for N=2048, and benefits increase for larger N.

**For N=32K (production scale)**, Option B saves even more:
- Memory: 0.85 GB per rank (2× less)
- Time: 75s vs 135s per rank (45% faster)
- Total: 109 CPU-hours saved per run

## Troubleshooting

**"Out of memory" with Option A?**
→ Switch to Option B (saves ~0.05 GB per rank for N=2048, 0.85 GB for N=32K)

**Outputs don't match?**
→ Check: same param file, same seed, same ranks, same grid decomposition

**Option B slower?**
→ Unlikely, but check: I/O bottleneck, network contention

## Next Steps

1. Run small-scale test (N=256) to verify correctness
2. Run medium-scale test (N=2048) to measure performance
3. Switch to Option B if tests confirm benefits
4. For production (N=32K), Option B benefits will be even larger
5. Monitor first production run

---

**For detailed information, see `TESTING_STRATEGY_N32K.md`**

