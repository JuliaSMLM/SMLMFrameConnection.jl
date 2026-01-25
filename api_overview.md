# SMLMFrameConnection API Overview

Frame-connection for 2D SMLM data: combines repeated localizations of blinking fluorophores into higher-precision localizations.

## Main Functions

### frameconnect
```julia
frameconnect(smld::BasicSMLD{T,Emitter2DFit{T}};
    nnearestclusters::Int=2,
    nsigmadev::Float64=5.0,
    maxframegap::Int=5,
    nmaxnn::Int=2
) where T -> (smld_connected, smld_preclustered, smld_combined, params)
```
Main entry point. Connects repeated localizations and combines them.

**Parameters:**
- `nnearestclusters`: Nearest preclusters for local density estimation
- `nsigmadev`: Sigma multiplier defining preclustering distance threshold (higher = larger connection radius)
- `maxframegap`: Maximum frame gap for temporal adjacency in preclusters
- `nmaxnn`: Maximum nearest-neighbors inspected for precluster membership

**Returns:**
- `smld_connected`: Input with `track_id` populated (localizations uncombined)
- `smld_preclustered`: Intermediate preclustering result (for debugging)
- `smld_combined`: **Main output** - combined high-precision localizations
- `params::ParamStruct`: Algorithm parameters (input + estimated photophysics)

### combinelocalizations
```julia
combinelocalizations(smld::BasicSMLD{T,Emitter2DFit{T}}) where T -> BasicSMLD
```
Combines emitters sharing the same `track_id` using MLE weighted mean. Use when `track_id` is already populated.

### defineidealFC
```julia
defineidealFC(smld::BasicSMLD{T,Emitter2DFit{T}};
    maxframegap::Int=5
) where T -> (smld_connected, smld_combined)
```
For simulated data where `track_id` indicates ground-truth emitter ID. Useful for validation/benchmarking.

## Types

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
- `dataset`: Dataset identifier (default: 1)
- `track_id`: Set to 0 for input (populated by algorithm)

## Output

Output `smld_combined` contains `Emitter2DFit` emitters with:
- Combined position via MLE weighted mean: `x = Σ(x/σ²) / Σ(1/σ²)`
- Reduced uncertainties: `σ = √(1/Σ(1/σ²))`
- Summed photons
- `track_id` indicating connection group

## Dependencies

- SMLMData.jl (v0.5+): BasicSMLD, Emitter2DFit types
- Hungarian.jl (v0.6-0.7): Linear assignment problem solver
- NearestNeighbors.jl: Spatial clustering
- Optim.jl: Parameter estimation

## Algorithm Reference

Schodt & Lidke (2021), "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Front. Bioinform.
https://doi.org/10.3389/fbinf.2021.724325
