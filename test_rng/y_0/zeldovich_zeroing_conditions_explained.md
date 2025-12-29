# Explanation of Zeldovich Zeroing Conditions

## Overview

The zeldovich.cpp code has three conditions that can zero out D values (lines 360-364). Understanding these conditions is crucial for matching the behavior in the hermitian code.

## The Three Conditions

```cpp
if ( (abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax)
     || (!param.CornerModes && k2>=k2_cutoff)
     || (param.qonemode && !(kx==param.one_mode[0] && ky==param.one_mode[1] && kz==param.one_mode[2])) ) {
    D = 0.0;
    if (ver == 2) nskip++;
}
```

---

## Condition 1: Nyquist Frequency Zeroing

```cpp
abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax
```

### Purpose
Zeros out modes at the Nyquist frequency (or k_cutoff-adjusted Nyquist) along any axis.

### Details
- `kmax = (ppdhalf * ik_cutoff + 0.5)` where `ik_cutoff = 1.0 / k_cutoff`
- For `k_cutoff = 1.0`: `kmax = 8` (Nyquist frequency for N=16)
- For `k_cutoff = 2.0`: `kmax = 4` (half Nyquist)
- **Why zero these?** Nyquist modes are self-conjugate and don't have separate partners. They must be zero to maintain proper Hermitian symmetry.

### Examples (N=16, k_cutoff=1.0)
- `k = (8, 0, 0)`: `abs(kx)=8 == kmax` → **Zeroed**
- `k = (0, 8, 0)`: `abs(ky)=8 == kmax` → **Zeroed**
- `k = (0, 0, 8)`: `abs(kz)=8 == kmax` → **Zeroed**
- `k = (6, 0, 6)`: All `abs(k) < 8` → **Not zeroed by this condition**

---

## Condition 2: k_cutoff Filtering (The Missing Condition!)

```cpp
!param.CornerModes && k2>=k2_cutoff
```

### Purpose
**Filters out high-wavenumber modes** to enable oversampling/convergence tests. This is the condition that's **missing in the hermitian code**.

### Details
- `k2_cutoff = (nyquist²) / (k_cutoff²) = (N/2)² / (k_cutoff²)`
- For N=16, k_cutoff=1.0: `k2_cutoff = 8² / 1.0² = 64`
- **When `CornerModes = 0` (default)**: All modes with `k² >= k2_cutoff` are zeroed
- **When `CornerModes = 1`**: These modes are NOT zeroed (allows "corner modes" with k > k_Nyquist)

### Why This Exists
From the README:
> "This is useful for doing convergence tests, e.g. run once with `PPD=64` and `ZD_k_cutoff = 1`, and again with `PPD=128` and `ZD_k_cutoff = 2`. This will produce two boxes with the exact same modes (although the PLT corrections will be slightly different), but the second box's modes are oversampled by a factor of two."

**Purpose**: Generate ICs with different resolutions but **identical modes** for convergence testing.

### Examples (N=16, k_cutoff=1.0, CornerModes=0)
- `k = (6, 0, 6)`: `k² = 72 >= 64` → **Zeroed** ✓ (This is why (6,6) differs!)
- `k = (7, 0, 4)`: `k² = 65 >= 64` → **Zeroed** ✓
- `k = (4, 0, 7)`: `k² = 65 >= 64` → **Zeroed** ✓
- `k = (5, 0, 5)`: `k² = 50 < 64` → **Not zeroed**
- `k = (8, 0, 0)`: `k² = 64 >= 64` → **Zeroed** (also caught by Condition 1)

### What Are CornerModes?
From `parameters.h`:
```cpp
int CornerModes;  // fill modes k > k_Ny. Default: 0.
```

- **CornerModes = 0 (default)**: Zero out all modes with `k² >= k2_cutoff`
- **CornerModes = 1**: Allow modes with `k > k_Nyquist` to be filled (corner modes in k-space)

**When would you use CornerModes=1?**
- Rarely used in practice
- Allows power beyond the Nyquist sphere to be generated
- Typically not desired for standard IC generation

---

## Condition 3: One-Mode Filtering (Debugging/Testing)

```cpp
param.qonemode && !(kx==param.one_mode[0] && ky==param.one_mode[1] && kz==param.one_mode[2])
```

### Purpose
**Debugging/testing feature**: Generate ICs with only ONE specific mode, zeroing all others.

### Details
- **When `qonemode = 0` (default)**: This condition is never true → No effect
- **When `qonemode > 0`**: Only the mode specified in `one_mode[3]` is kept, all others are zeroed

### Example
If `qonemode = 1` and `one_mode = {4, 0, 0}`:
- `k = (4, 0, 0)`: Matches `one_mode` → **Not zeroed** (kept)
- `k = (5, 0, 0)`: Doesn't match → **Zeroed**
- `k = (0, 0, 0)`: Doesn't match → **Zeroed**

### Use Case
From the README:
> "This is useful for automatically iterating through a series of wavevectors, for examining isotropy or Nyquist effects, for example."

**Purpose**: Test individual modes in isolation to verify code behavior.

---

## Summary Table

| Condition | Default | Purpose | When Active |
|-----------|---------|---------|-------------|
| **Nyquist** | Always | Zero self-conjugate modes | `abs(kx\|ky\|kz) == kmax` |
| **k_cutoff** | Active | Filter high-k modes | `CornerModes=0 && k² >= k2_cutoff` |
| **qonemode** | Inactive | Debug: keep one mode only | `qonemode > 0 && k != one_mode` |

---

## Impact on Current Differences

### The Missing Condition

The hermitian code is **missing Condition 2** (k_cutoff filtering). This causes differences for all coordinates where:
- `k² >= k2_cutoff = 64` (for N=16, k_cutoff=1.0)
- `CornerModes = 0` (default)

### Examples of Affected Coordinates

All 8 non-mirrored Y=0 differences have `k² >= 64`:
- (4,7): `k² = 65` → Zeroed in zeldovich, non-zero in hermitian
- (5,7): `k² = 74` → Zeroed in zeldovich, non-zero in hermitian
- (6,6): `k² = 72` → Zeroed in zeldovich, non-zero in hermitian
- (6,7): `k² = 85` → Zeroed in zeldovich, non-zero in hermitian
- (7,4): `k² = 65` → Zeroed in zeldovich, non-zero in hermitian
- (7,5): `k² = 74` → Zeroed in zeldovich, non-zero in hermitian
- (7,6): `k² = 85` → Zeroed in zeldovich, non-zero in hermitian
- (7,7): `k² = 98` → Zeroed in zeldovich, non-zero in hermitian

### Why This Matters

1. **Convergence Testing**: The k_cutoff feature enables comparing simulations at different resolutions with identical modes
2. **Physical Correctness**: High-k modes beyond k_cutoff should be zeroed to match zeldovich behavior
3. **RNG Consistency**: The `nskip++` ensures RNG alignment is maintained even when modes are zeroed

---

## Implementation Notes

### Required Changes in Hermitian Code

1. **Read k_cutoff parameter** (default: 1.0)
2. **Read CornerModes parameter** (default: 0)
3. **Calculate k2_cutoff**: `k2_cutoff = (Nhalf²) / (k_cutoff²)`
4. **Add zeroing condition**:
   ```c
   else if (!params->CornerModes && k2_int >= k2_cutoff) {
       D[0] = D[1] = 0.0;
       nskip++;  // Maintain RNG alignment
   }
   ```

### Parameter File

The parameter file should support:
- `ZD_k_cutoff` (default: 1.0)
- `ZD_CornerModes` (default: 0)

These are already in the zeldovich parameter file format, so the hermitian code should read them the same way.

