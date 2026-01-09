# SMLMFrameConnection

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl)

## Overview

SMLMFrameConnection performs frame-connection on localization microscopy data organized in an `SMLMData.BasicSMLD` structure with `Emitter2DFit` emitters (see [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl)).

Frame-connection combines repeated localizations of a single blinking event into a single higher-precision localization. This is done using the algorithm presented in [Schodt & Lidke 2021](https://doi.org/10.3389/fbinf.2021.724325). Connected localizations share the same `track_id` value in the output emitters.

## Installation

```julia
using Pkg
Pkg.add("SMLMFrameConnection")
```

## Basic Usage

```julia
using SMLMData, SMLMFrameConnection

# Create or load your BasicSMLD with Emitter2DFit emitters
# Emitters must have valid x, y, σ_x, σ_y, and frame fields

smld_connected, smld_preclustered, smld_combined, params = frameconnect(smld)
```

### Outputs

- `smld_connected`: Input with `track_id` updated to associate connected localizations
- `smld_preclustered`: Intermediate result with precluster associations
- `smld_combined`: **Main output** - combined higher-precision localizations
- `params`: Parameters used in the algorithm (user-defined and estimated)

### Parameters

```julia
smld_connected, smld_preclustered, smld_combined, params = frameconnect(smld;
    nnearestclusters = 2,   # Nearest preclusters for density estimation
    nsigmadev = 5.0,        # Sigma multiplier for preclustering threshold
    maxframegap = 5,        # Maximum frame gap for temporal adjacency
    nmaxnn = 2              # Maximum nearest-neighbors for precluster membership
)
```

## Input Requirements

The input `BasicSMLD` must contain `Emitter2DFit` emitters with:
- `x`, `y`: Position coordinates (microns)
- `σ_x`, `σ_y`: Position uncertainties (microns) - **required for MLE combination**
- `frame`: Frame number
- `dataset`: Dataset identifier (for multi-dataset processing)
- `photons`, `bg`, `σ_photons`, `σ_bg`: Photometry (combined in output)

## Example

```julia
using SMLMData, SMLMFrameConnection

# Create camera and emitters
cam = IdealCamera(1:512, 1:512, 0.1)  # 0.1 μm pixels
emitters = [
    Emitter2DFit{Float64}(5.0, 5.0, 1000.0, 10.0, 0.02, 0.02, 50.0, 2.0; frame=1),
    Emitter2DFit{Float64}(5.01, 5.01, 1200.0, 12.0, 0.02, 0.02, 60.0, 2.0; frame=2),
    Emitter2DFit{Float64}(5.02, 4.99, 1100.0, 11.0, 0.02, 0.02, 55.0, 2.0; frame=3),
]
smld = BasicSMLD(emitters, cam, 3, 1)

# Run frame connection
_, _, smld_combined, _ = frameconnect(smld; maxframegap=5)

# smld_combined now contains 1 emitter with improved precision
```

## Citation

David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Frontiers in Bioinformatics, 2021

https://doi.org/10.3389/fbinf.2021.724325
