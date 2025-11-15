# This file defines some struct types used in the FrameConnection package.

"""
    ParamStruct

Immutable structure of parameters needed for frame-connection.

# Fields
- `initialdensity`: Density of emitters at the start of the experiment (emitters/pixel^2).
                    See `estimatedensities()`.
- `nnearestclusters`: Number of nearest preclusters used for local density
                      estimates (default = 2). See `estimatedensities()`.
- `k_on`: Rate at which dark emitters convert to the visible state (1/frame).
          See `estimateparams()`.
- `k_off`: Rate at which visible emitters are converted to the reversible dark
           state (1/frame). See `estimateparams()`.
- `k_bleach`: Rate at which visible emitters are irreversibly photobleached (1/frame).
              See `estimateparams()`.
- `p_miss`: Probability of missing a localization of a visible emitter.
- `nsigmadev`: Multiplier of localization errors that defines a pre-clustering
               distance threshold (default = 5, in units of pixels). See `precluster()`.
- `maxframegap`: Maximum frame gap between temporally adjacent localizations in
                 a precluster (default = 5 frames). See `precluster()`.
- `nmaxnn`: Maximum number of nearest-neighbors inspected for precluster
            membership (default = 2). Ideally infinite, but not feasible for large data.
            See `precluster()`.
"""
struct ParamStruct
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

# Default constructor
ParamStruct() = ParamStruct(Float64[], 2, 0.0, 0.0, 0.0, 0.0, 5.0, 5, 2)