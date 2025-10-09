# InSAR DEM Toolbox

A modular, high-resolution InSAR processing pipeline for generating digital elevation models (DEMs) from dual-polarization single-look complex (SLC) data.

---

## 📁 Folder Structure

```
InSAR/
├── InSAR_DEM_Pipeline.m            # Main processing script
├── setup_insar_toolbox.m           # Toolbox path initializer
├── README.md                       # You are here
│
├── colormaps/                      # Custom colormaps (e.g., RdYlBu.csv)
└── InSAR_DEM_Modules/              # Modular namespaced functions
    ├── +insar/                     # Interferogram, unwrapping, DEM
    ├── +io/                        # I/O and file reading
    └── +utils/                     # Filtering, geometry, projection
```

---

## 🚀 Getting Started

1. Open MATLAB
2. Navigate to the `InSAR/` root directory
3. Run setup and pipeline:

```matlab
setup_insar_toolbox()
run InSAR_DEM_Pipeline
```

---

## 📦 Requirements

### MATLAB Version
- **R2020b or later** (R2021a+ preferred for Mapping Toolbox updates)

### Required Toolboxes
- **Image Processing Toolbox**
- **Mapping Toolbox**

### File Dependencies
Ensure the following are available in your dataset:
- Geocoded `.mli_geo.tif` files (for reference geometry)
- Binary `.slc` files (single-look complex SAR data)
- LiDAR-derived `.tif` DEM (for phase-to-height conversion)
- Flight trajectory files with azimuth time metadata

---

## 🧠 Authors

**Tate Meehan**  
Glaciology & Remote Sensing | [@WCP](https://www.wildernesscollectivepreserve.org)

**ChatGPT (code structuring & refactoring assistant)**  
OpenAI, July 2025

---

## 🛠 Contributions & Acknowledgments

This toolbox structure was inspired by scientific modular design best practices and built for operational use during the 2024–2025 WCP radar campaigns.

---

## 🔄 License

This project is currently internal/research-use only. Licensing terms to be defined.
