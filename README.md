# SMLMFrameConnection

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl)

Frame-connection for single-molecule localization microscopy: linking localizations from the same fluorophore blinking event across consecutive frames into single, higher-precision localizations. Uses spatiotemporal LAP assignment to optimally connect temporally adjacent detections based on spatial proximity and estimated blinking kinetics.

## Installation

```julia
using Pkg
Pkg.add("SMLMFrameConnection")
```

## Quick Start

```julia
using SMLMFrameConnection

# Frame connection on localization data
(combined, info) = frameconnect(smld)

# Output: combined high-precision localizations
println("$(info.n_input) → $(info.n_combined) localizations")
```

For complete SMLM workflows (detection + fitting + frame-connection + rendering), see [SMLMAnalysis.jl](https://github.com/JuliaSMLM/SMLMAnalysis.jl).

## Configuration

`frameconnect()` accepts keyword arguments or a `FrameConnectConfig` struct:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_density_neighbors` | 2 | Nearest preclusters for local density estimation |
| `max_sigma_dist` | 5.0 | Sigma multiplier for preclustering distance threshold |
| `max_frame_gap` | 5 | Maximum frame gap for temporal adjacency |
| `max_neighbors` | 2 | Maximum nearest-neighbors for precluster membership |

```julia
# Keyword form (most common)
(combined, info) = frameconnect(smld; max_frame_gap=10, max_sigma_dist=3.0)

# Config struct form (reusable settings)
config = FrameConnectConfig(max_frame_gap=10, max_sigma_dist=3.0)
(combined, info) = frameconnect(smld, config)
```

**Parameter guidance:** Default values work well for standard dSTORM/PALM data. For dense samples, reduce `max_sigma_dist` to 3.0. For long dark states (dSTORM), increase `max_frame_gap` to 10-20.

## Output Format

`frameconnect()` returns `(combined::BasicSMLD, info::FrameConnectInfo)`.

| Output | Description |
|--------|-------------|
| `combined` | High-precision combined localizations (main output) |
| `info.connected` | Input with `track_id` assigned (per-frame data with labels) |
| `info.n_input` | Number of input localizations |
| `info.n_tracks` | Number of tracks formed |
| `info.n_combined` | Number of output localizations |
| `info.k_on` | Estimated on rate (1/frame) |
| `info.k_off` | Estimated off rate (1/frame) |
| `info.k_bleach` | Estimated bleach rate (1/frame) |
| `info.p_miss` | Probability of missed detection |
| `info.initial_density` | Density estimate per cluster (emitters/μm²) |
| `info.elapsed_s` | Wall time (seconds) |
| `info.algorithm` | Algorithm used (`:lap`) |
| `info.n_preclusters` | Number of preclusters formed |

### Combination Method

Connected localizations are combined using maximum likelihood estimation (MLE) weighted mean:
- Position: `x_combined = Σ(x/σ²) / Σ(1/σ²)` (inverse-variance weighted)
- Uncertainty: `σ_combined = √(1/Σ(1/σ²))` ≈ `σ_individual / √n`
- Photons: summed across connected localizations

## Algorithm Pipeline

1. **Precluster**: Spatiotemporal clustering using KDTree nearest-neighbor search within `max_sigma_dist * σ` distance and `max_frame_gap` frames
2. **Estimate parameters**: Fit photophysics rates (k_on, k_off, k_bleach, p_miss) from precluster statistics using Optim.jl NelderMead
3. **Connect via LAP**: Build cost matrix from spatial separation and photophysics likelihoods, solve with Hungarian.jl
4. **Combine**: MLE weighted mean using full 2×2 covariance (precision-weighted)

## Example

```julia
using SMLMFrameConnection

# Create camera (pixel_ranges, pixel_size in μm)
cam = IdealCamera(1:512, 1:512, 0.1)

# Create emitters representing the same molecule blinking across 3 frames
emitters = [
    Emitter2DFit{Float64}(
        5.00, 5.00,    # x, y position (μm)
        1000.0, 10.0,  # photons, background
        0.02, 0.02, 0.0,  # σ_x, σ_y, σ_xy
        50.0, 2.0,     # σ_photons, σ_bg
        1, 1, 0, 1     # frame, dataset, track_id, id
    ),
    Emitter2DFit{Float64}(5.01, 5.01, 1200.0, 12.0, 0.02, 0.02, 0.0, 60.0, 2.0, 2, 1, 0, 2),
    Emitter2DFit{Float64}(5.02, 4.99, 1100.0, 11.0, 0.02, 0.02, 0.0, 55.0, 2.0, 3, 1, 0, 3),
]

smld = BasicSMLD(emitters, cam, 3, 1)

# Run frame connection
(combined, info) = frameconnect(smld)

# Result: Combined uncertainty ≈ σ_individual / √n_connected
println("$(info.n_input) → $(info.n_combined) localizations in $(info.elapsed_s)s")
```

## Utility Functions

### combinelocalizations
```julia
smld_combined = combinelocalizations(smld)
```
Combines emitters with the same `track_id`. Use when you have pre-labeled data.

### defineidealFC
```julia
smld_connected, smld_combined = defineidealFC(smld; max_frame_gap=5)
```
For simulated data where `track_id` contains ground-truth emitter IDs. Validates frame-connection performance against known truth.

## Algorithm Reference

> Schodt, D.J. and Lidke, K.A. "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem." *Frontiers in Bioinformatics*, 2021. [DOI: 10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)

## Related Packages

- **[SMLMAnalysis.jl](https://github.com/JuliaSMLM/SMLMAnalysis.jl)** - Complete SMLM workflow (detection + fitting + frame-connection + rendering)
- **[SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl)** - Core data types for SMLM
- **[GaussMLE.jl](https://github.com/JuliaSMLM/GaussMLE.jl)** - GPU-accelerated Gaussian PSF fitting
- **[SMLMSim.jl](https://github.com/JuliaSMLM/SMLMSim.jl)** - SMLM data simulation

## License

MIT License - see [LICENSE](LICENSE) file for details.
