# Options A vs B Testing

This directory contains all scripts and outputs for testing and comparing Options A and B for particle IC writing.

## Directory Structure

```
tests/options_ab/
├── README.md                    # This file
├── test_option_a_N256.pbs      # PBS script for Option A (N=256)
├── test_option_b_N256.pbs      # PBS script for Option B (N=256)
├── test_option_a_N2048.pbs     # PBS script for Option A (N=2048)
├── test_option_b_N2048.pbs     # PBS script for Option B (N=2048)
├── submit_all_tests.sh         # Submit all tests at once
├── compare_option_outputs.sh    # Compare outputs from both options
├── test_options_A_vs_B_N32K.sh # Alternative bash-based test script
├── output_option_a/            # Output directory for Option A (created by jobs)
├── output_option_b/            # Output directory for Option B (created by jobs)
├── test_option_*_N256.out      # PBS output files (N=256)
├── test_option_*_N256.log      # Execution logs (N=256)
├── test_option_*_N2048.out     # PBS output files (N=2048)
└── test_option_*_N2048.log      # Execution logs (N=2048)
```

## Quick Start

### 1. Submit Tests

From the project root:
```bash
cd tests/options_ab
./submit_all_tests.sh --small-only    # Submit N=256 tests only
./submit_all_tests.sh --medium-only   # Submit N=2048 tests only
./submit_all_tests.sh                 # Submit all tests
```

Or from project root:
```bash
cd tests/options_ab && ./submit_all_tests.sh --small-only
```

### 2. Monitor Jobs

```bash
qstat -u $USER
```

### 3. Compare Results (After Jobs Complete)

**Compare outputs (file counts, sizes, data):**
```bash
cd tests/options_ab
./compare_option_outputs.sh output_option_a output_option_b

# Or with data verification
./compare_option_outputs.sh output_option_a output_option_b --verify-data
```

**Compare timing (performance):**
```bash
cd tests/options_ab
./compare_timing.sh 256    # Compare N=256 timing
./compare_timing.sh 2048   # Compare N=2048 timing
./compare_timing.sh         # Compare all available tests
```

## Output Locations

- **Option A outputs**: `tests/options_ab/output_option_a/`
- **Option B outputs**: `tests/options_ab/output_option_b/`
- **PBS output files**: `tests/options_ab/test_option_*.out`
- **Execution logs**: `tests/options_ab/test_option_*.log`

## Notes

- All PBS scripts change to the project root before running
- Output directories are created relative to the project root
- Log files are saved in this directory for easy access

