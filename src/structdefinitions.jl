# This file defines some struct types used in the FrameConnection package.

"""
    CalibrationConfig

Configuration for optional uncertainty calibration within frame connection.

Calibration analyzes frame-to-frame jitter within connected tracks to estimate
a motion variance (`σ_motion²`) and CRLB scale factor (`k²`), then applies corrected
uncertainties before track combination: `Σ_corrected = σ_motion² I + k² Σ_CRLB`.

# Fields
- `clamp_k_to_one::Bool=true`: Clamp k_scale to minimum of 1.0 (CRLB is a lower bound)
- `filter_high_chi2::Bool=false`: Filter tracks with high chi² pairs before fitting
- `chi2_filter_threshold::Float64=6.0`: Chi² threshold for track filtering (per pair)
"""
Base.@kwdef struct CalibrationConfig
    clamp_k_to_one::Bool = true
    filter_high_chi2::Bool = false
    chi2_filter_threshold::Float64 = 6.0
end

"""
    CalibrationResult

Output diagnostics from uncertainty calibration. Returned in `FrameConnectInfo.calibration`
when calibration is enabled.

The calibration model fits: `observed_var = A + B * CRLB_var` where
`A = σ_motion²` (additive motion/vibration variance) and `B = k²` (CRLB scale factor).

# Fields
- `sigma_motion_nm::Float64`: Estimated motion std dev in nm (√A × 1000)
- `k_scale::Float64`: CRLB scale factor (√B, or clamped √B if clamp_k_to_one)
- `A::Float64`: Fit intercept (σ_motion² in μm²)
- `B::Float64`: Fit slope (k²)
- `A_sigma::Float64`: Standard error of A
- `B_sigma::Float64`: Standard error of B
- `r_squared::Float64`: R² of the WLS fit
- `mean_chi2::Float64`: Mean chi² across all frame-to-frame pairs
- `n_pairs::Int`: Number of frame-to-frame pairs used in fit
- `n_tracks_used::Int`: Number of tracks contributing pairs
- `n_tracks_filtered::Int`: Number of tracks removed by chi² filter
- `bin_centers::Vector{Float64}`: Bin centers (CRLB variance) for diagnostic plots
- `bin_observed::Vector{Float64}`: Bin observed variance for diagnostic plots
- `frame_shifts::Dict{Int, Vector{NTuple{2,Float64}}}`: Per-dataset (dx,dy) shifts for jitter plots
- `calibration_applied::Bool`: Whether calibration was actually applied (false on fallback)
- `warning::String`: Warning message if calibration fell back (empty if OK)
"""
struct CalibrationResult
    sigma_motion_nm::Float64
    k_scale::Float64
    A::Float64
    B::Float64
    A_sigma::Float64
    B_sigma::Float64
    r_squared::Float64
    mean_chi2::Float64
    n_pairs::Int
    n_tracks_used::Int
    n_tracks_filtered::Int
    bin_centers::Vector{Float64}
    bin_observed::Vector{Float64}
    frame_shifts::Dict{Int, Vector{NTuple{2,Float64}}}
    calibration_applied::Bool
    warning::String
end

"""
    FrameConnectConfig

Configuration parameters for frame connection algorithm.

# Fields
- `n_density_neighbors::Int=2`: Number of nearest preclusters used for local density
                             estimates (see `estimatedensities`)
- `max_sigma_dist::Float64=5.0`: Multiplier of localization errors that defines a
                            pre-clustering distance threshold (see `precluster`)
- `max_frame_gap::Int=5`: Maximum frame gap between temporally adjacent localizations
                        in a precluster (see `precluster`)
- `max_neighbors::Int=2`: Maximum number of nearest-neighbors inspected for precluster
                   membership (see `precluster`)
- `calibration::Union{CalibrationConfig, Nothing}=nothing`: Optional uncertainty calibration

# Example
```julia
# Without calibration (default)
(combined, info) = frameconnect(smld)

# With calibration
config = FrameConnectConfig(calibration=CalibrationConfig())
(combined, info) = frameconnect(smld, config)

# With calibration + chi2 filtering
config = FrameConnectConfig(
    calibration=CalibrationConfig(filter_high_chi2=true, chi2_filter_threshold=4.0)
)
(combined, info) = frameconnect(smld, config)
```
"""
Base.@kwdef struct FrameConnectConfig <: AbstractSMLMConfig
    n_density_neighbors::Int = 2
    max_sigma_dist::Float64 = 5.0
    max_frame_gap::Int = 5
    max_neighbors::Int = 2
    calibration::Union{CalibrationConfig, Nothing} = nothing
end

"""
    ParamStruct(initial_density::Vector{Float64},
        k_on::Float64, k_off::Float64, k_bleach::Float64, p_miss::Float64,
        max_sigma_dist::Float64, max_frame_gap::Int, n_density_neighbors::Int)

Structure of parameters needed for frame-connection.

# Fields
- `initial_density`: Density of emitters at the start of the experiment.
                    (see estimatedensities()) (emitters/μm²)
- `n_density_neighbors`: Number of nearest preclusters used for local density 
                      estimates. (default = 2)(see estimatedensities())
- `k_on`: Rate at which dark emitters convert to the visible state (1/frame).
- `k_off`: Rate at which visible emitters convert to the reversible dark state (1/frame).
- `k_bleach`: Rate at which visible emitters are irreversibly photobleached (1/frame).
- `p_miss`: Probability of missing a localization of a visible emitter.
- `max_sigma_dist`: Multiplier of localization errors that defines a pre-clustering
               distance threshold. (default = 5)(see precluster())(unitless)
- `max_frame_gap`: Maximum frame gap between temporally adjacent localizations in
                 a precluster. (default = 5)(see precluster())(frames)
- `max_neighbors`: Maximum number of nearest-neighbors inspected for precluster
            membership.  Ideally, this would be set to inf, but that's not
            feasible for most data. (default = 2)(see precluster())

Note: k_on, k_off, and k_bleach are transition rates, not duty cycle fractions.
The duty cycle (fraction of time in ON state) is k_on/(k_on + k_off). For typical
dSTORM, k_on << k_off gives low duty cycle (mostly dark, brief blinks).
"""
mutable struct ParamStruct
    initial_density::Vector{Float64}
    n_density_neighbors::Int
    k_on::Float64
    k_off::Float64
    k_bleach::Float64
    p_miss::Float64
    max_sigma_dist::Float64
    max_frame_gap::Int
    max_neighbors::Int
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