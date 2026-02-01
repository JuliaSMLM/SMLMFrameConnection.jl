```@meta
CurrentModule = SMLMFrameConnection
```

# SMLMFrameConnection

Frame-connection for 2D single molecule localization microscopy (SMLM) data. Combines repeated localizations of blinking fluorophores into higher-precision localizations using the algorithm from [Schodt & Lidke 2021](https://doi.org/10.3389/fbinf.2021.724325).

## Installation

```julia
using Pkg
Pkg.add("SMLMFrameConnection")
```

## Quick Start

```julia
using SMLMData, SMLMFrameConnection

# Run frame connection on your BasicSMLD with Emitter2DFit emitters
(combined, info) = frameconnect(smld)

# combined is the main output - higher precision localizations
# info contains connected SMLD, statistics, and estimated photophysics
```

## Input Requirements

Input `BasicSMLD` must contain `Emitter2DFit` emitters with:

**Required:**
- `x`, `y`: Position coordinates (microns)
- `σ_x`, `σ_y`: Position uncertainties (microns) - must be > 0 for MLE weighting
- `frame`: Frame number (1-based)

**Optional:**
- `photons`, `bg`: Photometry (summed in output)
- `σ_xy`: Position covariance (propagated in output)
- `dataset`: Dataset identifier (default: 1)

## Algorithm Overview

1. **Preclustering**: Groups spatially and temporally adjacent localizations into candidate clusters
2. **Parameter Estimation**: Estimates fluorophore blinking kinetics (`k_on`, `k_off`, `k_bleach`) and emitter density from the data
3. **Frame Connection**: Uses Linear Assignment Problem (LAP) to optimally assign localizations to emitters based on spatial proximity and estimated photophysics
4. **Combination**: Combines connected localizations using MLE weighted mean with full covariance propagation

## Outputs

```julia
(combined, info) = frameconnect(smld)
```

| Output | Type | Description |
|--------|------|-------------|
| `combined` | `BasicSMLD` | **Main output** - combined high-precision localizations |
| `info` | `ConnectInfo` | Metadata: connected SMLD, statistics, photophysics |

### ConnectInfo Fields

| Field | Description |
|-------|-------------|
| `info.connected` | Original localizations with `track_id` labels |
| `info.n_input` | Number of input localizations |
| `info.n_tracks` | Number of tracks formed |
| `info.n_combined` | Number of output localizations |
| `info.n_preclusters` | Number of preclusters |
| `info.k_on`, `k_off`, `k_bleach` | Estimated photophysics |
| `info.p_miss` | Estimated miss probability |
| `info.elapsed_ns` | Processing time (nanoseconds) |
| `info.algorithm` | Algorithm used (`:lap`) |

## Parameters

```julia
frameconnect(smld;
    nnearestclusters = 2,   # Clusters for density estimation
    nsigmadev = 5.0,        # Distance threshold multiplier
    maxframegap = 5,        # Max frame gap for connections
    nmaxnn = 2              # Nearest neighbors for preclustering
)
```

- `nsigmadev`: Higher values allow connections over larger distances
- `maxframegap`: Increase for dyes with long dark states (dSTORM: 10-20)

## Example

```julia
using SMLMData, SMLMFrameConnection

# Run frame connection
(combined, info) = frameconnect(smld; maxframegap=10)

# Access results
println("Combined $(info.n_input) → $(info.n_combined) localizations")
println("Formed $(info.n_tracks) tracks")
println("Estimated k_on=$(info.k_on), k_off=$(info.k_off)")

# Access connected (uncombined) localizations with track_id
connected_smld = info.connected
```

## API Reference

```@index
```

```@autodocs
Modules = [SMLMFrameConnection]
```

## Citation

David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Frontiers in Bioinformatics, 2021. [DOI: 10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)
