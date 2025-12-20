# Why F(k) = delta(k) * i * k_x / k²

## The Physics: Zeldovich Approximation

In the Zeldovich approximation, the displacement field **q(x)** is related to the density field **δ(x)** through:

1. **Poisson's equation** relates gravitational potential to density
2. **Displacement** is the gradient of the potential

## Step-by-Step Derivation

### Step 1: Poisson's Equation in Real Space

```
∇²Φ = (3/2) * Ω_m * H₀² * δ
```

Where:
- Φ = gravitational potential
- δ = density perturbation
- Ω_m = matter density parameter
- H₀ = Hubble constant

### Step 2: Fourier Transform

In Fourier space, the Laplacian becomes multiplication by -k²:

```
∇² → -k²
```

So Poisson's equation in Fourier space:
```
-k² * Φ(k) = (3/2) * Ω_m * H₀² * δ(k)
```

Solving for potential:
```
Φ(k) = -(3/2) * Ω_m * H₀² * δ(k) / k²
```

### Step 3: Displacement from Potential

The displacement field is the **gradient** of the potential:

```
q(x) = -∇Φ(x)
```

### Step 4: Gradient in Fourier Space

**Key insight:** Taking a gradient in real space corresponds to **multiplication by i*k** in Fourier space!

```
∇ → i*k  (in Fourier space)
```

So:
```
q(k) = -∇Φ(k) = -i*k * Φ(k)
```

### Step 5: Combine Everything

Substitute Φ(k) from Step 2:

```
q(k) = -i*k * [-(3/2) * Ω_m * H₀² * δ(k) / k²]
     = i*k * (3/2) * Ω_m * H₀² * δ(k) / k²
     = i * (k/k²) * (3/2) * Ω_m * H₀² * δ(k)
```

### Step 6: Component-Wise

For the x-component:
```
q_x(k) = i * (k_x/k²) * (3/2) * Ω_m * H₀² * δ(k)
```

The normalization factor (3/2) * Ω_m * H₀² is typically absorbed into the power spectrum normalization, so we get:

```
q_x(k) = i * (k_x/k²) * δ(k)
```

Or, in Zeldovich's notation:
```
F(k) = i * (k_x/k²) * D(k)
```

Where D(k) = δ(k) is the density Fourier coefficient.

## Why the "i" (Imaginary Unit)?

**The "i" comes from the gradient operator in Fourier space!**

- Real space gradient: `∇f(x)` 
- Fourier space: `i*k * f(k)`

This is a fundamental property of Fourier transforms:
```
F[∇f(x)] = i*k * F[f(x)]
```

The imaginary unit appears because:
- Taking a derivative multiplies by i*k in Fourier space
- This ensures the phase relationships are correct

## Why k_x / k²?

- **k_x**: The x-component of the wavevector (direction of displacement)
- **k²**: Normalization factor from Poisson's equation (1/k² comes from the Laplacian)

The ratio k_x/k² gives:
- **Direction**: k_x points in the x-direction
- **Magnitude**: 1/k² provides the correct scaling (longer wavelengths have larger displacements)

## Visual Intuition

For a density perturbation at wavevector k:
- The displacement is **perpendicular to the density wave** (in the direction of k)
- The magnitude scales as 1/k² (longer wavelengths → larger displacements)
- The phase relationship (the "i") ensures the displacement leads/lags the density correctly

## Summary

```
F(k) = i * (k_x/k²) * δ(k)
```

- **i**: From gradient operator in Fourier space
- **k_x**: X-component of wavevector (direction)
- **k²**: Normalization from Poisson's equation
- **δ(k)**: Density Fourier coefficient
