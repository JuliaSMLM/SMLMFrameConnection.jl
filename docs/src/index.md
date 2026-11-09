```@meta
CurrentModule = SMLMFrameConnection
```

# SMLMFrameConnection

Frame-connection for 2D single molecule localization microscopy (SMLM) data: linking localizations from the same fluorophore blinking event across consecutive frames into single, higher-precision localizations. Uses spatiotemporal LAP assignment to optimally connect temporally adjacent detections based on spatial proximity and estimated blinking kinetics. See [Schodt & Lidke 2021](https://doi.org/10.3389/fbinf.2021.724325).

## Installation

```julia
using Pkg
Pkg.add("SMLMFrameConnection")
```

## Quick Start

```julia
using SMLMFrameConnection

# Run frame connection on your BasicSMLD with Emitter2DFit emitters
(combined, info) = frameconnect(smld)

# combined is the main output - higher precision localizations
# info contains track assignments and algorithm metadata
```

## Input Requirements

Input `BasicSMLD` must contain `Emitter2DFit` emitters with:

**Required:**
- `x`, `y`: Position coordinates (microns), finite
- `σ_x`, `σ_y`: Position uncertainties (microns), finite and > 0
- `σ_xy`: Position covariance (μm²); the covariance must be positive definite, `|σ_xy| < σ_x σ_y`
- `frame`: Frame number; for `frameconnect`, at least 1, and `smld.n_frames` is at least the largest frame, or 0 (unknown)

**Optional:**
- `photons`, `bg`: Photometry (summed in output)
- `dataset`: Dataset identifier (default: 1)

`frameconnect` refuses input that breaks these with an `ArgumentError` naming the
localization and field; `combinelocalizations` checks the same except the frame rules (it
only groups by track and takes each track's earliest frame). Localizations that carry `z`
(3D) are refused the same way by both; frame connection is 2D only.

Localizations may come in any order (for example SMLMSim's output, grouped by emitter):
`frameconnect` works on a copy sorted by every field, so the result is the same for any
input order, and `info.connected` keeps the caller's order.

`combinelocalizations` combines localizations that share a `track_id` within a dataset.
`track_id` 0 means unlinked: each such localization passes through on its own, as does a
track of one localization.

## Algorithm Overview

1. **Preclustering**: Groups spatially and temporally adjacent localizations into candidate clusters
2. **Parameter Estimation**: Estimates fluorophore blinking kinetics (`k_on`, `k_off`, `k_bleach`) and emitter density from the data
3. **Frame Connection**: Uses Linear Assignment Problem (LAP) to optimally assign localizations to emitters based on spatial proximity and estimated photophysics
4. **Combination**: Combines connected localizations using MLE weighted mean for improved precision

## Outputs

| Output | Description |
|--------|-------------|
| `combined` | **Main output** - combined high-precision localizations |
| `info.connected` | Original localizations with `track_id` assigned |
| `info.n_tracks` | Number of tracks formed |
| `info.elapsed_s` | Wall time in seconds |

## Parameters

```julia
frameconnect(smld;
    n_density_neighbors = 2,   # Clusters for density estimation
    max_sigma_dist = 5.0,        # Distance threshold multiplier
    max_frame_gap = 5,        # Max frame gap for connections
    max_neighbors = 2,              # Nearest neighbors for preclustering
    track_length = nothing      # Inclusive (min,max) locs-per-track filter; nothing disables
)
```

- `max_sigma_dist`: Higher values allow connections over larger distances
- `max_frame_gap`: Increase for dyes with long dark states (dSTORM: 10-20)
- `track_length`: `(2.0, Inf)` drops single-frame blinks; `(2.0, 50.0)` also drops over-long tracks (fiducials, sticky strands). Count dropped is reported in `info.n_filtered`

## API Reference

```@index
```

```@autodocs
Modules = [SMLMFrameConnection]
```

## Citation

David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Frontiers in Bioinformatics, 2021. [DOI: 10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)
