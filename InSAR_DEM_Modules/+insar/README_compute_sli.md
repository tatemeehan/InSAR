
# compute_sli.m — Calibrated Single Look Image with Optional Radiometric Corrections

This MATLAB function computes amplitude, power, and decibel-scaled single-look images from complex SLC data, with options for:
- Corner reflector-based σ⁰ calibration
- Pixel area normalization
- Radiometric terrain corrections: σ⁰ (incidence), γ⁰ (slope), β⁰ (range)

## 📥 Inputs

```matlab
[amp, pow, db] = compute_sli(slc, P, a, lambda, pixelArea, incAngle, slopeDeg, mode)
```

- `slc`        : Complex single look complex (SLC) image (MxN)
- `P`          : Measured return power at corner reflector (linear scale)
- `a`          : Corner reflector side length (default = 1 m)
- `lambda`     : Radar wavelength (m) (default = 0.3 / 1.3 for L-band)
- `pixelArea`  : Ground pixel area in m² (scalar or MxN raster)
- `incAngle`   : Incidence angle in degrees (MxN raster)
- `slopeDeg`   : Terrain slope in degrees (MxN raster)
- `mode`       : Radiometric correction type: `'sigma0'`, `'gamma0'`, `'beta0'` (default: `'sigma0'`)

## 📤 Outputs

- `amp`        : Amplitude image (|SLC|)
- `pow`        : Calibrated power image (linear scale)
- `db`         : Calibrated power image in decibels

## 📐 Radiometric Corrections

- **σ⁰ (sigma0)**: Normalized radar cross section (backscatter per unit area perpendicular to radar beam)
- **γ⁰ (gamma0)**: Normalized to ground plane (removes slope distortions)
- **β⁰ (beta0)** : Normalized to slant range pixel (no terrain correction)

Requires:
- Incidence angle raster for σ⁰/γ⁰
- Terrain slope raster for γ⁰
- Pixel area (optional, but used if provided)

## 📘 Usage Example

```matlab
[amp, pow, db] = compute_sli(slcData, P_CR, 1, 0.23, 1, incRaster, slopeRaster, 'gamma0');
```

## 🧑‍🔬 Author
ChatGPT + Tate Meehan  
Last updated: July 2025
