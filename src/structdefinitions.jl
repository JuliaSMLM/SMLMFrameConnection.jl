# This file defines some struct types used in the FrameConnection package.

"""
    FrameConnectConfig

Configuration parameters for frame connection algorithm.

# Fields
- `nnearestclusters::Int=2`: Number of nearest preclusters used for local density
                             estimates (see `estimatedensities`)
- `nsigmadev::Float64=5.0`: Multiplier of localization errors that defines a
                            pre-clustering distance threshold (see `precluster`)
- `maxframegap::Int=5`: Maximum frame gap between temporally adjacent localizations
                        in a precluster (see `precluster`)
- `nmaxnn::Int=2`: Maximum number of nearest-neighbors inspected for precluster
                   membership (see `precluster`)

# Example
```julia
# Using default config
config = FrameConnectConfig()
(combined, info) = frameconnect(smld, config)

# Custom config
config = FrameConnectConfig(maxframegap=10, nsigmadev=3.0)
(combined, info) = frameconnect(smld, config)

# Kwargs form (equivalent to Config form)
(combined, info) = frameconnect(smld; maxframegap=10, nsigmadev=3.0)
```
"""
Base.@kwdef struct FrameConnectConfig <: AbstractSMLMConfig
    nnearestclusters::Int = 2
    nsigmadev::Float64 = 5.0
    maxframegap::Int = 5
    nmaxnn::Int = 2
end

"""
    ParamStruct(initialdensity::Vector{Float64},
        k_on::Float64, k_off::Float64, k_bleach::Float64, p_miss::Float64,
        nsigmadev::Float64, maxframegap::Int, nnearestclusters::Int)

Structure of parameters needed for frame-connection.

# Fields
- `initialdensity`: Density of emitters at the start of the experiment.
                    (see estimatedensities()) (emitters/μm²)
- `nnearestclusters`: Number of nearest preclusters used for local density 
                      estimates. (default = 2)(see estimatedensities())
- `k_on`: Rate at which dark emitters convert to the visible state (1/frame).
- `k_off`: Rate at which visible emitters convert to the reversible dark state (1/frame).
- `k_bleach`: Rate at which visible emitters are irreversibly photobleached (1/frame).
- `p_miss`: Probability of missing a localization of a visible emitter.
- `nsigmadev`: Multiplier of localization errors that defines a pre-clustering
               distance threshold. (default = 5)(see precluster())(unitless)
- `maxframegap`: Maximum frame gap between temporally adjacent localizations in
                 a precluster. (default = 5)(see precluster())(frames)
- `nmaxnn`: Maximum number of nearest-neighbors inspected for precluster
            membership.  Ideally, this would be set to inf, but that's not
            feasible for most data. (default = 2)(see precluster())

Note: k_on, k_off, and k_bleach are transition rates, not duty cycle fractions.
The duty cycle (fraction of time in ON state) is k_on/(k_on + k_off). For typical
dSTORM, k_on << k_off gives low duty cycle (mostly dark, brief blinks).
"""
mutable struct ParamStruct
    initialdensity::Vector{Float64}
    nnearestclusters::Int
    k_on::Float64
    k_off::Float64
    k_bleach::Float64
    p_miss::Float64
    nsigmadev::Float64
    maxframegap::Int
    nmaxnn::Int
end
ParamStruct() = ParamStruct([], 2, 0.0, 0.0, 0.0, 0.0, 5.0, 5, 2)

"""
    to_emitter2dfit(e::AbstractEmitter, track_id::Int) -> Emitter2DFit

Convert any AbstractEmitter to Emitter2DFit with specified track_id.
Used internally to support different input emitter types while producing
consistent Emitter2DFit output.
"""
function to_emitter2dfit(e::SMLMData.AbstractEmitter, track_id::Int)
    T = typeof(e.x)
    SMLMData.Emitter2DFit{T}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
        e.frame, e.dataset, track_id, e.id
    )
end