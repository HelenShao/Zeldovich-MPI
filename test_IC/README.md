# Brute Force Trace Comparison (N=8)

This directory contains scripts and workflows for brute force tracing of matrix values through the computation pipeline, comparing `hermitian_3d_matrix_production` and `zeldovich-PLT`.

## Purpose

Trace and compare matrix values at two critical points:
1. **Before FFT**: After generation in Fourier space (D, F, G, H arrays)
2. **After FFT**: After 3D FFT in real space (final displacement/velocity fields)

## Directory Structure

```
test_IC/
├── param_N8_trace.par          # Parameter file for N=8 trace test
├── workflow_trace.txt           # Step-by-step workflow
├── README.md                    # This file
├── run_zeldovich_trace.sh       # Script to run zeldovich-PLT with trace output
├── extract_matrix_dumps.py      # Extract matrix values from logs
├── compare_matrices.py          # Compare matrix dumps between codes
├── generate_trace_report.py     # Generate summary report
├── write_ic_trace.sh            # Reassemble .bin files and write particle ICs
├── hermitian_output/            # Output from hermitian_3d_matrix_production
│   ├── hermitian_N8_trace.log
│   ├── matrix_before_fft.txt
│   └── matrix_after_fft.txt
├── zeldovich_output/            # Output from zeldovich-PLT
│   ├── zeldovich_N8_trace.log
│   ├── matrix_before_fft.txt
│   └── matrix_after_fft.txt
├── matrix_dumps/                # Raw matrix dumps (if needed)
├── comparisons/                 # Comparison results
│   ├── before_fft_comparison.txt
│   ├── after_fft_comparison.txt
│   └── trace_report.txt
└── particle_ics/                 # Particle IC files (for end-to-end verification)
```

## Quick Start

1. Follow the workflow in `workflow_trace.txt`
2. Check comparison results in `comparisons/`
3. Review summary report: `comparisons/trace_report.txt`

## Matrix Format

Matrices are dumped in the following format:
- **Before FFT**: Fourier space, organized by (Y, X, Z) with array index
  - Format: `Y=<y> X=<x> Z=<z> Array=<array_idx> Re=<real> Im=<imag>`
- **After FFT**: Real space, organized by (Z, Y, X) with array index
  - Format: `Z=<z> Y=<y> X=<x> Array=<array_idx> Re=<real> Im=<imag>`

## Notes

- N=8 is chosen for manageable matrix size (8³ = 512 elements per array)
- Both codes use the same random seed (ZD_Seed = 4) for reproducibility
- PLT is disabled for direct comparison
- k_cutoff = 1.0, CornerModes = 0 (same as k2_cutoff test)

