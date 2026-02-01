# SMLMFrameConnection API Overview

Frame-connection for 2D SMLM data: combines repeated localizations of blinking fluorophores into higher-precision localizations.

## Exports

4 exports: `frameconnect`, `combinelocalizations`, `defineidealFC`, `ConnectInfo`, `ParamStruct`

## Main Functions

### frameconnect
```julia
(combined, info) = frameconnect(smld::BasicSMLD{T,Emitter2DFit{T}};
    nnearestclusters::Int=2,
    nsigmadev::Float64=5.0,
    maxframegap::Int=5,
    nmaxnn::Int=2
) where T
```
Main entry point. Connects repeated localizations and combines them.

**Parameters:**
- `nnearestclusters`: Nearest preclusters for local density estimation
- `nsigmadev`: Sigma multiplier defining preclustering distance threshold (higher = larger connection radius)
- `maxframegap`: Maximum frame gap for temporal adjacency in preclusters
- `nmaxnn`: Maximum nearest-neighbors inspected for precluster membership

**Returns tuple:**
- `combined::BasicSMLD`: **Main output** - combined high-precision localizations
- `info::ConnectInfo`: Metadata including connected SMLD, statistics, and estimated photophysics

### combinelocalizations
```julia
combined = combinelocalizations(smld::BasicSMLD{T,Emitter2DFit{T}}) where T -> BasicSMLD
```
Combines emitters sharing the same `track_id` using MLE weighted mean with full covariance propagation. Use when `track_id` is already populated.

### defineidealFC
```julia
(connected, combined) = defineidealFC(smld::BasicSMLD{T,Emitter2DFit{T}};
    maxframegap::Int=5
) where T
```
For simulated data where `track_id` indicates ground-truth emitter ID. Useful for validation/benchmarking.

## Types

### ConnectInfo{T}
```julia
struct ConnectInfo{T}
    connected::BasicSMLD{T}       # Input with track_id assigned (uncombined)
    n_input::Int                  # Number of input localizations
    n_tracks::Int                 # Number of tracks formed
    n_combined::Int               # Number of output localizations
    k_on::Float64                 # Estimated on rate (1/frame)
    k_off::Float64                # Estimated off rate (1/frame)
    k_bleach::Float64             # Estimated bleach rate (1/frame)
    p_miss::Float64               # Probability of missed detection
    initialdensity::Vector{Float64}  # Density estimate per cluster (emitters/μm²)
    elapsed_ns::UInt64            # Processing time in nanoseconds
    algorithm::Symbol             # Algorithm used (:lap)
    n_preclusters::Int            # Number of preclusters formed
end
```

### ParamStruct
```julia
mutable struct ParamStruct
    initialdensity::Vector{Float64}  # Emitter density per cluster (emitters/μm²)
    nnearestclusters::Int            # Clusters for density estimation
    k_on::Float64                    # Rate: dark → visible (1/frame)
    k_off::Float64                   # Rate: visible → dark (1/frame)
    k_bleach::Float64                # Rate: photobleaching (1/frame)
    p_miss::Float64                  # Probability of missing localization when on
    nsigmadev::Float64               # Preclustering threshold multiplier
    maxframegap::Int                 # Max frame gap in preclusters
    nmaxnn::Int                      # Max nearest-neighbors for preclustering
end
```

**Photophysics terms:**
- `k_on`: Transition rate from dark (non-fluorescent) to visible state
- `k_off`: Transition rate from visible to dark state
- `k_bleach`: Irreversible photobleaching rate
- `p_miss`: Probability localization algorithm fails to detect an "on" fluorophore

## Input Requirements

Input `BasicSMLD` must contain `Emitter2DFit` emitters.

**Required fields:**
- `x`, `y`: Position (microns)
- `σ_x`, `σ_y`: Position uncertainties (microns) - must be > 0 for MLE
- `frame`: Frame number (1-based)

**Optional fields:**
- `photons`, `σ_photons`: Photon count (summed in output)
- `bg`, `σ_bg`: Background
- `σ_xy`: Position covariance (microns², propagated in output)
- `dataset`: Dataset identifier (default: 1)
- `track_id`: Set to 0 for input (populated by algorithm)

## Output

Output `combined` contains `Emitter2DFit` emitters with:
- Combined position via MLE weighted mean with full 2x2 covariance
- Properly propagated uncertainties including σ_xy
- Summed photons
- `track_id` indicating connection group

## Example

```julia
using SMLMData, SMLMFrameConnection

# Load or create SMLD
smld = ...

# Run frame connection
(combined, info) = frameconnect(smld; maxframegap=10)

# Access results
println("Combined $(info.n_input) → $(info.n_combined) localizations")
println("Estimated k_on=$(info.k_on), k_off=$(info.k_off)")
println("Processing time: $(info.elapsed_ns / 1e9) seconds")

# Access connected (uncombined) localizations
connected_smld = info.connected
```

## Dependencies

- SMLMData.jl (v0.6+): BasicSMLD, Emitter2DFit types
- Hungarian.jl (v0.6-0.7): Linear assignment problem solver
- NearestNeighbors.jl: Spatial clustering
- Optim.jl: Parameter estimation

## Algorithm Reference

Schodt & Lidke (2021), "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Front. Bioinform.
https://doi.org/10.3389/fbinf.2021.724325
