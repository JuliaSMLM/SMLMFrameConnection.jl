# SMLMFrameConnection

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl)

## Overview

SMLMFrameConnection performs **frame-connection** on 2D localization microscopy data: combining repeated localizations of a single blinking fluorophore into a single higher-precision localization.

Uses the spatiotemporal clustering algorithm from [Schodt & Lidke 2021](https://doi.org/10.3389/fbinf.2021.724325).

## Installation

```julia
using Pkg
Pkg.add("SMLMFrameConnection")
```

## Quick Start

```julia
using SMLMData, SMLMFrameConnection

# Run frame connection on your data
smld_connected, smld_preclustered, smld_combined, params = frameconnect(smld)

# smld_combined is the main output - higher precision localizations
```

## Input Requirements

Input must be a `BasicSMLD{T, Emitter2DFit{T}}` from [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl).

**Required fields** (algorithm fails without these):
- `x`, `y`: Position coordinates in microns
- `σ_x`, `σ_y`: Position uncertainties in microns (used for MLE weighting; must be > 0)
- `frame`: Frame number (1-based integer)

**Optional fields** (combined in output if present):
- `photons`, `σ_photons`: Photon count and uncertainty (summed across connected localizations)
- `bg`, `σ_bg`: Background and uncertainty
- `dataset`: Dataset identifier (defaults to 1; for multi-ROI or multi-acquisition data)

## Complete Example

```julia
using SMLMData, SMLMFrameConnection

# Create camera (pixel_ranges, pixel_size in μm)
cam = IdealCamera(1:512, 1:512, 0.1)

# Create emitters representing the same molecule blinking across 3 frames
# Constructor: Emitter2DFit{T}(x, y, photons, bg, σ_x, σ_y, σ_photons, σ_bg; frame, dataset, track_id)
emitters = [
    Emitter2DFit{Float64}(
        5.0, 5.0,      # x, y position (μm)
        1000.0, 10.0,  # photons, background
        0.02, 0.02,    # σ_x, σ_y uncertainties (μm)
        50.0, 2.0;     # σ_photons, σ_bg
        frame=1
    ),
    Emitter2DFit{Float64}(5.01, 5.01, 1200.0, 12.0, 0.02, 0.02, 60.0, 2.0; frame=2),
    Emitter2DFit{Float64}(5.02, 4.99, 1100.0, 11.0, 0.02, 0.02, 55.0, 2.0; frame=3),
]

# Create SMLD: BasicSMLD(emitters, camera, n_frames, n_datasets)
smld = BasicSMLD(emitters, cam, 3, 1)

# Run frame connection
_, _, smld_combined, _ = frameconnect(smld)

# Result: localizations connected based on spatial/temporal proximity
# Combined uncertainty: σ_combined ≈ σ_individual / √n_connected
```

## Outputs Explained

```julia
smld_connected, smld_preclustered, smld_combined, params = frameconnect(smld)
```

| Output | Description | When to use |
|--------|-------------|-------------|
| `smld_combined` | **Main output.** Combined high-precision localizations | Standard analysis |
| `smld_connected` | Original localizations with `track_id` populated | When you need per-frame data with connection labels |
| `smld_preclustered` | Intermediate preclustering result | Debugging, algorithm tuning |
| `params` | Estimated photophysics + input parameters | Inspecting estimated k_on, k_off, density |

## Parameters

```julia
frameconnect(smld;
    nnearestclusters = 2,   # Clusters used for local density estimation
    nsigmadev = 5.0,        # Distance threshold = nsigmadev × localization uncertainty
    maxframegap = 5,        # Max frames between connected localizations
    nmaxnn = 2              # Nearest neighbors checked during preclustering
)
```

**Parameter guidance:**
- `nsigmadev`: Higher values allow connections over larger distances. Default (5.0) works for typical SMLM data. Reduce for dense samples.
- `maxframegap`: Set based on expected blinking duration. For dSTORM with long dark states, increase to 10-20.
- Defaults work well for standard dSTORM/PALM data with typical blinking kinetics.

## Estimated Parameters (ParamStruct)

The algorithm estimates fluorophore photophysics from your data:

| Field | Description |
|-------|-------------|
| `k_on` | Rate of transitioning from dark to visible state (1/frame) |
| `k_off` | Rate of transitioning from visible to dark state (1/frame) |
| `k_bleach` | Photobleaching rate (1/frame) |
| `p_miss` | Probability of missing a localization when fluorophore is on |
| `initialdensity` | Estimated emitter density per cluster (emitters/μm²) |

## Combination Method

Connected localizations are combined using **maximum likelihood estimation (MLE) weighted mean**:
- Position: inverse-variance weighted average → `x_combined = Σ(x/σ²) / Σ(1/σ²)`
- Uncertainty: `σ_combined = √(1/Σ(1/σ²))` ≈ `σ_individual / √n`
- Photons: summed across connected localizations

## Utility Functions

### combinelocalizations
```julia
smld_combined = combinelocalizations(smld)
```
Combines emitters that share the same `track_id`. Use when you have pre-labeled data.

### defineidealFC
```julia
smld_connected, smld_combined = defineidealFC(smld; maxframegap=5)
```
For **simulated data** where `track_id` already contains ground-truth emitter IDs. Useful for validating frame-connection performance against known truth.

## Citation

David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Frontiers in Bioinformatics, 2021

https://doi.org/10.3389/fbinf.2021.724325
