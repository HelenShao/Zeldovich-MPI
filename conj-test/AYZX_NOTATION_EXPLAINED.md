# AYZX Notation and Array Structure in Zeldovich-PLT

## AYZX Macro Definition

```cpp
#define AYZX(_slab, _a, _y, _z, _x) \
    _slab[(int64_t) (_x) + array.ppd * ((_z) + array.ppd * ((_a) + array.narray * (_y)))]
```

**Parameters:**
- `_slab`: Pointer to the data slab
- `_a`: Array index (0, 1, 2, 3)
- `_y`: Y coordinate
- `_z`: Z coordinate  
- `_x`: X coordinate

**Memory Layout:** `[Y][Array][Z][X]`
- X is stride-1 (fastest varying)
- Z is next
- Array index is next
- Y is stride-1 for FFT along Y direction

**Formula:** `x + ppd * (z + ppd * (array_idx + narray * y))`

## What Each Array Contains

### Array 0: Density + X-Displacement
```cpp
AYZX(slab, 0, yres, z, x) = D + I * F;
```

**In Fourier Space:**
- `D` = density field (complex)
- `F` = x-displacement field (complex) = `i * (k_x/k²) * D`

**After 3D IFFT (Real Space):**
- `real(Array[0])` = density field (real-valued)
- `imag(Array[0])` = x-displacement field (real-valued)

### Array 1: Y-Displacement + Z-Displacement
```cpp
AYZX(slab, 1, yres, z, x) = G + I * H;
```

**In Fourier Space:**
- `G` = y-displacement field (complex) = `i * (k_y/k²) * D`
- `H` = z-displacement field (complex) = `i * (k_z/k²) * D`

**After 3D IFFT (Real Space):**
- `real(Array[1])` = y-displacement field (real-valued)
- `imag(Array[1])` = z-displacement field (real-valued)

### Array 2: X-Velocity (if PLT enabled)
```cpp
AYZX(slab, 2, yres, z, x) = 0. + I * F * f;
```
- `f` = PLT growth rate
- `real(Array[2])` = 0
- `imag(Array[2])` = x-velocity = `F * f`

### Array 3: Y-Velocity + Z-Velocity (if PLT enabled)
```cpp
AYZX(slab, 3, yres, z, x) = G * f + I * H * f;
```
- `real(Array[3])` = y-velocity = `G * f`
- `imag(Array[3])` = z-velocity = `H * f`

## Field Relationships

All displacement fields are derived from density:
```
F = i * (k_x/k²) * D  // X-displacement
G = i * (k_y/k²) * D  // Y-displacement
H = i * (k_z/k²) * D  // Z-displacement
```

This comes from:
1. Poisson's equation: `∇²Φ = density`
2. Displacement: `q = -∇Φ`
3. In Fourier space: `q(k) = i * (k/k²) * density(k)`

## How Fields Are Extracted (from output.cpp)

```cpp
// Density
dens = real(YX(slab1, y, x));  // From Array[0].real

// Displacements
pos[0] = imag(YX(slab1, y, x)) * norm;  // X from Array[0].imag
pos[1] = real(YX(slab2, y, x)) * norm;  // Y from Array[1].real
pos[2] = imag(YX(slab2, y, x)) * norm;  // Z from Array[1].imag
```

Where:
- `slab1` = Array[0] (density + x-displacement)
- `slab2` = Array[1] (y-displacement + z-displacement)

## Summary

| Array | Real Part | Imaginary Part |
|-------|-----------|----------------|
| 0     | Density   | X-displacement |
| 1     | Y-displacement | Z-displacement |
| 2     | 0         | X-velocity (if PLT) |
| 3     | Y-velocity | Z-velocity (if PLT) |

**Key Point:** Each array stores TWO real-valued fields packed into one complex array structure. After 3D IFFT, both components are real numbers (not complex), but they're stored in the real and imaginary slots of the complex array.

