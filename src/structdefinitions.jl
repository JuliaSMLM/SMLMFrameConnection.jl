# This file defines some struct types used in the FrameConnection package.

"""
    CalibrationConfig

Configuration for optional uncertainty calibration within frame connection.

Calibration analyzes frame-to-frame jitter within connected tracks to estimate
a motion variance (`σ_motion²`) and CRLB scale factor (`k²`), then applies corrected
uncertainties before track combination: `Σ_corrected = σ_motion² I + k² Σ_CRLB`.

# Fields
- `clamp_k_to_one::Bool=true`: Hold k_scale ≥ 1.0 (CRLB is a lower bound); when the fit gives
  k < 1, the motion term is refit with k = 1
- `filter_high_chi2::Bool=false`: Filter tracks with high chi² pairs before fitting. The chi² here
  uses the uncalibrated CRLB, so this filter removes the tracks that show a too-small CRLB and
  pulls `k_scale` toward 1; the fit itself already leaves out mismatched pairs (see
  [`CalibrationResult`](@ref)).
- `chi2_filter_threshold::Float64=6.0`: Chi² threshold for track filtering (per pair)
"""
Base.@kwdef struct CalibrationConfig <: SMLMData.AbstractSMLMConfig
    clamp_k_to_one::Bool = true
    filter_high_chi2::Bool = false
    chi2_filter_threshold::Float64 = 6.0
end

"""
    CalibrationResult

Output diagnostics from uncertainty calibration. Returned in `FrameConnectInfo.calibration`
when calibration is enabled.

The calibration model fits: `observed_var = A + B * CRLB_var` where
`observed_var` is the per-axis variance of the difference of two consecutive
localizations, `CRLB_var` the sum of their per-axis CRLB variances, so the intercept is
`A = 2σ_motion²` (each localization carries its own motion/vibration jitter) and
`B = k²` (CRLB scale factor).

The fit leaves out pairs whose chi² under the current fit exceeds the chi²(2) 99.7% point
(mismatched pairs), corrects the truncated bin means for the variance the cut removes, and
iterates to convergence, so `k_scale` and `sigma_motion_nm` are unbiased when the model holds.

# Fields
- `sigma_motion_nm::Float64`: Applied per-localization motion std dev in nm (√(A_applied/2) × 1000)
- `k_scale::Float64`: Applied CRLB scale factor (√B_applied)
- `A::Float64`: Intercept of the free two-parameter fit (pair-difference variance 2σ_motion², in μm²); can be negative
- `B::Float64`: Slope of the free two-parameter fit (k²); below 1 when the clamp binds. The
  line `A + B * bin_centers` goes through `bin_observed`
- `A_sigma::Float64`: Standard error of A
- `B_sigma::Float64`: Standard error of B
- `r_squared::Float64`: R² of the WLS fit
- `mean_chi2::Float64`: Mean chi² across all frame-to-frame pairs
- `n_pairs::Int`: Number of frame-to-frame pairs used in the final fit (after the chi² cut)
- `n_tracks_used::Int`: Number of tracks contributing pairs
- `n_tracks_filtered::Int`: Number of tracks removed by chi² filter
- `bin_centers::Vector{Float64}`: Bin centers (CRLB variance) for diagnostic plots
- `bin_observed::Vector{Float64}`: Bin observed variance (truncation-corrected) of the free fit, for diagnostic plots
- `bin_counts::Vector{Int}`: Number of pairs per bin (parallel to bin_centers/bin_observed); use to grey/drop low-statistics bins in diagnostic plots
- `frame_shifts::Dict{Int, Vector{NTuple{2,Float64}}}`: Per-dataset (dx,dy) shifts for jitter plots
- `calibration_applied::Bool`: Whether calibration was actually applied (false on fallback)
- `warning::String`: Why calibration fell back (too few pairs or bins, or nothing
  determined), or, on an applied calibration, that motion and CRLB scale were not
  separable and which one-parameter fit was applied, or that the fit did not converge
  (empty if OK)
- `A_applied::Float64`: Intercept applied, in μm² (≥ 0): the free fit's, or a refit with a
  bound held (A ≥ 0; B ≥ 1 with `clamp_k_to_one`), or a one-parameter fit when the free
  fit cannot separate motion from the CRLB scale. 0 on fallback
- `B_applied::Float64`: Slope applied (k², > 0, ≥ 1 with the clamp). 1 on fallback
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
    bin_counts::Vector{Int}
    frame_shifts::Dict{Int, Vector{NTuple{2,Float64}}}
    calibration_applied::Bool
    warning::String
    A_applied::Float64
    B_applied::Float64
end

# The 17-field form (before A_applied and B_applied existed): the applied values follow
# from sigma_motion_nm and k_scale.
function CalibrationResult(sigma_motion_nm, k_scale, A, B, A_sigma, B_sigma, r_squared,
                           mean_chi2, n_pairs, n_tracks_used, n_tracks_filtered,
                           bin_centers, bin_observed, bin_counts, frame_shifts,
                           calibration_applied, warning)
    return CalibrationResult(sigma_motion_nm, k_scale, A, B, A_sigma, B_sigma, r_squared,
                             mean_chi2, n_pairs, n_tracks_used, n_tracks_filtered,
                             bin_centers, bin_observed, bin_counts, frame_shifts,
                             calibration_applied, warning,
                             2 * (sigma_motion_nm / 1000)^2, k_scale^2)
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
- `track_length::Union{Tuple{Float64,Float64}, Nothing}=nothing`: Inclusive `(min, max)`
                   range on the number of localizations a track must contain to be emitted
                   as a combined localization (`lo <= nperID <= hi`). Counts localizations
                   sharing a `track_id`, not temporal frame-span. `nothing` (default) disables
                   filtering. Use `(2.0, Inf)` to drop single-frame blinks, or `(2.0, 50.0)`
                   to also drop over-long tracks (e.g. fiducials, sticky docking strands in
                   dense DNA-PAINT). Mirrors the `(min, max)` filter idiom used elsewhere in
                   the ecosystem. See `combinelocalizations`.
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
    track_length::Union{Tuple{Float64,Float64}, Nothing} = nothing
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
"""
    _check_2d(smld, fname)

Throw an `ArgumentError` if any localization in `smld` carries a z coordinate. Frame
connection is 2D: it builds `Emitter2DFit` output and would silently drop z and σ_z.
"""
function _check_2d(smld::BasicSMLD, fname::AbstractString)
    if any(e -> hasproperty(e, :z), smld.emitters)
        throw(ArgumentError("$fname: the localizations carry z " *
            "($(eltype(smld.emitters))); frame connection is 2D only and would drop " *
            "z and σ_z. Pass 2D localizations (for example Emitter2DFit)."))
    end
    return nothing
end

"""
    _check_inputs(smld, fname; frames = true)

Throw an `ArgumentError` naming the localization and field when an input would hang
the LAP or give wrong numbers without an error: a non-finite position; σ_x or σ_y not
finite and > 0 (a zero σ gives NaN combined positions, a NaN σ NaN costs); a
covariance σ_xy that is not finite or leaves Σ not positive definite; and, with
`frames`, a frame below 1 (frames are 1-based) or an `n_frames` smaller than the largest
frame (`n_frames = 0` means unknown and is allowed).
"""
function _check_inputs(smld::BasicSMLD, fname::AbstractString; frames::Bool = true)
    maxframe = 0
    for (i, e) in enumerate(smld.emitters)
        (isfinite(e.x) && isfinite(e.y)) ||
            throw(ArgumentError("$fname: localization $i has a non-finite position " *
                                "(x = $(e.x), y = $(e.y))"))
        (isfinite(e.σ_x) && isfinite(e.σ_y) && e.σ_x > 0 && e.σ_y > 0) ||
            throw(ArgumentError("$fname: localization $i has σ_x = $(e.σ_x), " *
                                "σ_y = $(e.σ_y); both must be finite and > 0"))
        (isfinite(e.σ_xy) && e.σ_xy^2 < e.σ_x^2 * e.σ_y^2) ||
            throw(ArgumentError("$fname: localization $i has σ_xy = $(e.σ_xy); the " *
                                "covariance must satisfy σ_xy² < σ_x²·σ_y²"))
        if frames
            e.frame >= 1 ||
                throw(ArgumentError("$fname: localization $i has frame $(e.frame); " *
                                    "frames are 1-based"))
            maxframe = max(maxframe, e.frame)
        end
    end
    if frames && smld.n_frames > 0 && smld.n_frames < maxframe
        throw(ArgumentError("$fname: n_frames = $(smld.n_frames) is smaller than the " *
                            "largest frame ($maxframe)"))
    end
    return nothing
end
