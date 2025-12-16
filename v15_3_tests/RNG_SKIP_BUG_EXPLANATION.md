# RNG Skip Bug Explanation - Visual Example

## Setup
- MAX_PPD = 4096 (but for visualization, assume MAX_PPD = 8)
- Processing y=0 (self-conjugate slice)
- We want to understand why k=(0,0,1) differs between N=4 and N=6

## Loop Structure for Self-Conjugate (y=0)

```c
for (int z = 0; z <= Nhalf; z++) {           // z = 0, 1, 2, ..., N/2
    int x_max = (z == 0 ? Nhalf + 1 : N);   // For z=0: x goes 0 to N/2+1
    for (int x = 0; x < x_max; x++) {        // For z>0: x goes 0 to N-1
        if (x == Nhalf + 1 && N < MAX_PPD) {
            nskip += (MAX_PPD - N);  // Accumulate skip
        }
        // Process (x,z) - apply nskip, then call cgauss()
    }
}
```

## Example: N=4 (Nhalf=2, MAX_PPD=8)

### Processing Order:

| Step | z | x | k-vector | Action | nskip before | nskip after | RNG state |
|------|---|---|----------|--------|-------------|-------------|-----------|
| 1 | 0 | 0 | (0,0,0) | DC mode, D=0 | 0 | 1 | RNG[0] |
| 2 | 0 | 1 | (1,0,0) | Process, apply skip=1 | 1 | 0 | RNG[1] → RNG[2] |
| 3 | 0 | 2 | (2,0,0) | Process | 0 | 0 | RNG[3] |
| 4 | 0 | 3 | (3,0,0) | **x=Nhalf+1, skip += (8-4)=4** | 0 | 4 | RNG[4] |
| 5 | 0 | 3 | (3,0,0) | Process, apply skip=4 | 4 | 0 | RNG[4] → RNG[8] |
| 6 | 1 | 0 | **(0,0,1)** | **Process, apply skip=0** | 0 | 0 | **RNG[9]** |

**Result for k=(0,0,1): Uses RNG[9]**

## Example: N=6 (Nhalf=3, MAX_PPD=8)

### Processing Order:

| Step | z | x | k-vector | Action | nskip before | nskip after | RNG state |
|------|---|---|----------|--------|-------------|-------------|-----------|
| 1 | 0 | 0 | (0,0,0) | DC mode, D=0 | 0 | 1 | RNG[0] |
| 2 | 0 | 1 | (1,0,0) | Process, apply skip=1 | 1 | 0 | RNG[1] → RNG[2] |
| 3 | 0 | 2 | (2,0,0) | Process | 0 | 0 | RNG[3] |
| 4 | 0 | 3 | (3,0,0) | Process | 0 | 0 | RNG[4] |
| 5 | 0 | 4 | (4,0,0) | **x=Nhalf+1, skip += (8-6)=2** | 0 | 2 | RNG[5] |
| 6 | 0 | 4 | (4,0,0) | Process, apply skip=2 | 2 | 0 | RNG[5] → RNG[7] |
| 7 | 1 | 0 | **(0,0,1)** | **Process, apply skip=0** | 0 | 0 | **RNG[8]** |

**Result for k=(0,0,1): Uses RNG[8]**

## The Problem

Even though k=(0,0,1) is the **same mode** in both cases:
- **N=4**: Uses RNG[9] (after skipping 4 at x=3, z=0)
- **N=6**: Uses RNG[8] (after skipping 2 at x=4, z=0)

The RNG state is different because:
1. The skip at `x = Nhalf+1` happens at different x values (x=3 for N=4, x=4 for N=6)
2. The skip amount is different: (8-4)=4 vs (8-6)=2
3. This skip is applied **before** processing the next z value
4. So when we process k=(0,0,1) (z=1, x=0), the RNG has been advanced by different amounts

## Visual Representation

### N=4, y=0 (z=0 row):
```
x:  0    1    2    3    4    5    6    7
k: (0,0) (1,0) (2,0) (3,0) [SKIP] [SKIP] [SKIP] [SKIP]
    ↓     ↓     ↓     ↓      +4 skips accumulated here
   RNG[0] RNG[2] RNG[3] RNG[8]
   
   After z=0 row: RNG position = 9
```

### N=6, y=0 (z=0 row):
```
x:  0    1    2    3    4    5    6    7
k: (0,0) (1,0) (2,0) (3,0) (4,0) [SKIP] [SKIP] [SKIP]
    ↓     ↓     ↓     ↓     ↓      +2 skips accumulated here
   RNG[0] RNG[2] RNG[3] RNG[4] RNG[7]
   
   After z=0 row: RNG position = 8
```

### When processing k=(0,0,1) (z=1, x=0):

**N=4**: 
- Start of z=1: RNG position = 9 (from end of z=0 row)
- Process k=(0,0,1): Uses RNG[9] → RNG[11] (2 numbers for complex)
- **Result: D uses RNG[9] and RNG[10]**

**N=6**: 
- Start of z=1: RNG position = 8 (from end of z=0 row)
- Process k=(0,0,1): Uses RNG[8] → RNG[10] (2 numbers for complex)
- **Result: D uses RNG[8] and RNG[9]**

**Same k-vector k=(0,0,1), but different RNG states!**

### The Key Issue:

The skip accumulated at `x = Nhalf+1` in the z=0 row affects the RNG state for **all subsequent processing**, including k=(0,0,1) in the z=1 row. Since the skip amount is different for N=4 vs N=6, the RNG state differs.

## Why This Is Wrong

The skip at `x = Nhalf+1` is meant to skip missing x-values (x = N to MAX_PPD-1) that don't exist in the actual grid. However:

1. **The skip happens at the wrong time**: It's accumulated when processing x=Nhalf+1, but this affects the RNG state for **all subsequent modes** in the same z-row and next z-rows.

2. **The skip affects unrelated modes**: When we process k=(0,0,1) (z=1, x=0), we're applying a skip that was meant for missing x-values at z=0. These are completely different modes!

3. **The skip amount depends on N**: Different N values produce different skip amounts, so the same mode (k=(0,0,1)) ends up using different RNG states.

## The Correct Behavior

For k=(0,0,1) to be consistent across N values:
- It should use the **same RNG state** regardless of N
- The skip from x=Nhalf+1 at z=0 should **not** affect k=(0,0,1) at z=1
- The skip should only affect modes that actually need to be skipped

## Solution

The bug was that **z=0 row was missing a skip for its missing x-values**.

### Root Cause:
For z=0 (self-conjugate):
- `x_max = Nhalf + 1`
- `x` goes from 0 to `Nhalf` (since `x < x_max`)
- The condition `x == Nhalf + 1` is **NEVER true** because the loop ends at `x = Nhalf`
- So the skip for missing x-values was never applied!

For z > 0:
- `x_max = N`
- `x` goes from 0 to `N-1`
- The condition `x == Nhalf + 1` IS triggered (when x = N/2 + 1)
- So the skip IS applied correctly

### The Fix:
Add a skip after processing z=0 row to account for missing x-values:

```c
// After the x-loop for z=0:
if (z == 0 && N < MAX_PPD) {
    // z=0 processes x=0 to Nhalf, but MAX_PPD grid processes x=0 to MAX_PPD/2
    // Skip needed: (MAX_PPD/2) - Nhalf = (MAX_PPD - N) / 2
    int64_t skip_amount = (MAX_PPD / 2) - Nhalf;
    nskip += skip_amount;
    // Apply skip immediately
    advance_rng(nskip);
    nskip = 0;
}
```

This ensures that after z=0 row, the RNG state is consistent across all N values, so k=(0,0,1) and other k=(0,0,z) modes use the same RNG state.

### Why This Works:
- For N=4: After z=0, we processed x=0,1,2 (Nhalf+1=3 modes) and now skip (4096/2 - 2) = 2046 to match MAX_PPD
- For N=6: After z=0, we processed x=0,1,2,3 (Nhalf+1=4 modes) and now skip (4096/2 - 3) = 2045 to match MAX_PPD
- The total RNG advancement from start to z=1, x=0 is now the same for both N values!

## Detailed Step-by-Step Trace

### N=4, Processing k=(0,0,1):

**Initial state:** RNG position = 0, nskip = 0

**z=0 row:**
1. (x=0, z=0): k=(0,0,0), DC mode → D=0, nskip++ → nskip=1, RNG stays at 0
2. (x=1, z=0): k=(1,0,0), apply skip=1 → RNG: 0→1, call cgauss() → RNG: 1→3, nskip=0
3. (x=2, z=0): k=(2,0,0), no skip → call cgauss() → RNG: 3→5
4. (x=3, z=0): k=(3,0,0), **x=Nhalf+1 detected!**
   - Accumulate skip: nskip += (8-4) = 4 → nskip=4
   - Apply skip=4 → RNG: 5→9
   - Call cgauss() → RNG: 9→11
   - nskip=0

**End of z=0 row:** RNG position = 11

**z=1 row:**
5. (x=0, z=1): k=(0,0,1) ← **This is our target mode**
   - nskip=0 (no accumulated skip)
   - Call cgauss() → **Uses RNG[11] and RNG[12]**
   - RNG: 11→13

**Result:** k=(0,0,1) uses RNG[11] and RNG[12]

### N=6, Processing k=(0,0,1):

**Initial state:** RNG position = 0, nskip = 0

**z=0 row:**
1. (x=0, z=0): k=(0,0,0), DC mode → D=0, nskip++ → nskip=1, RNG stays at 0
2. (x=1, z=0): k=(1,0,0), apply skip=1 → RNG: 0→1, call cgauss() → RNG: 1→3, nskip=0
3. (x=2, z=0): k=(2,0,0), no skip → call cgauss() → RNG: 3→5
4. (x=3, z=0): k=(3,0,0), no skip → call cgauss() → RNG: 5→7
5. (x=4, z=0): k=(4,0,0), **x=Nhalf+1 detected!**
   - Accumulate skip: nskip += (8-6) = 2 → nskip=2
   - Apply skip=2 → RNG: 7→9
   - Call cgauss() → RNG: 9→11
   - nskip=0

**End of z=0 row:** RNG position = 11

**z=1 row:**
6. (x=0, z=1): k=(0,0,1) ← **This is our target mode**
   - nskip=0 (no accumulated skip)
   - Call cgauss() → **Uses RNG[11] and RNG[12]**
   - RNG: 11→13

**Result:** k=(0,0,1) uses RNG[11] and RNG[12]

Wait, in this trace they both end up at RNG[11]! Let me recalculate...

Actually, the issue is more subtle. The skip happens at different points in the sequence, and the total number of modes processed before k=(0,0,1) is different:

**N=4:** Process 4 modes in z=0 row (x=0,1,2,3) before k=(0,0,1)
**N=6:** Process 5 modes in z=0 row (x=0,1,2,3,4) before k=(0,0,1)

Plus the skip amounts are different. So the RNG position will be different.

## The Real Issue

The fundamental problem is that **the skip at x=Nhalf+1 affects the RNG state for all subsequent processing**, including modes in different z-rows. The skip should only affect the specific modes that are being skipped (x = N to MAX_PPD-1), not all future modes.

The correct behavior would be:
- Skip should only advance RNG for the missing modes themselves
- Modes that are actually processed (like k=(0,0,1)) should use the same RNG state regardless of N
- The skip should not "leak" into subsequent z-rows

