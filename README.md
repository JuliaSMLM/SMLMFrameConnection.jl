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
(combined, info) = frameconnect(smld)

# combined is the main output - higher precision localizations
# info contains connected SMLD, statistics, and estimated photophysics
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
# Constructor: Emitter2DFit{T}(x, y, photons, bg, σ_x, σ_y, σ_xy, σ_photons, σ_bg, frame, dataset, track_id, id)
emitters = [
    Emitter2DFit{Float64}(5.0, 5.0, 1000.0, 10.0, 0.02, 0.02, 0.0, 50.0, 2.0, 1, 1, 0, 1),
    Emitter2DFit{Float64}(5.01, 5.01, 1200.0, 12.0, 0.02, 0.02, 0.0, 60.0, 2.0, 2, 1, 0, 2),
    Emitter2DFit{Float64}(5.02, 4.99, 1100.0, 11.0, 0.02, 0.02, 0.0, 55.0, 2.0, 3, 1, 0, 3),
]

# Create SMLD: BasicSMLD(emitters, camera, n_frames, n_datasets)
smld = BasicSMLD(emitters, cam, 3, 1)

# Run frame connection
(combined, info) = frameconnect(smld)

# Result: localizations connected based on spatial/temporal proximity
println("Combined $(info.n_input) localizations into $(info.n_combined)")
println("Formed $(info.n_tracks) tracks from $(info.n_preclusters) preclusters")

# Combined uncertainty: σ_combined ≈ σ_individual / √n_connected
```

## Outputs Explained

```julia
(combined, info) = frameconnect(smld)
```

| Output | Type | Description |
|--------|------|-------------|
| `combined` | `BasicSMLD` | **Main output.** Combined high-precision localizations |
| `info` | `ConnectInfo` | Metadata: connected SMLD, statistics, estimated photophysics |

### ConnectInfo Fields

| Field | Type | Description |
|-------|------|-------------|
| `info.connected` | `BasicSMLD` | Original localizations with `track_id` populated |
| `info.n_input` | `Int` | Number of input localizations |
| `info.n_tracks` | `Int` | Number of tracks formed |
| `info.n_combined` | `Int` | Number of output localizations |
| `info.n_preclusters` | `Int` | Number of preclusters |
| `info.k_on` | `Float64` | Estimated on rate (1/frame) |
| `info.k_off` | `Float64` | Estimated off rate (1/frame) |
| `info.k_bleach` | `Float64` | Estimated bleach rate (1/frame) |
| `info.p_miss` | `Float64` | Estimated miss probability |
| `info.initialdensity` | `Vector{Float64}` | Density estimate per cluster (emitters/μm²) |
| `info.elapsed_ns` | `UInt64` | Processing time in nanoseconds |
| `info.algorithm` | `Symbol` | Algorithm used (`:lap`) |

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

## Combination Method

Connected localizations are combined using **maximum likelihood estimation (MLE) weighted mean** with full covariance propagation:
- Position: inverse-variance weighted average using 2x2 covariance matrix
- Uncertainty: properly propagated including σ_xy correlation
- Photons: summed across connected localizations

## Utility Functions

### combinelocalizations
```julia
combined = combinelocalizations(smld)
```
Combines emitters that share the same `track_id`. Use when you have pre-labeled data.

### defineidealFC
```julia
(connected, combined) = defineidealFC(smld; maxframegap=5)
```
For **simulated data** where `track_id` already contains ground-truth emitter IDs. Useful for validating frame-connection performance against known truth.

## Citation

David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Frontiers in Bioinformatics, 2021

https://doi.org/10.3389/fbinf.2021.724325
