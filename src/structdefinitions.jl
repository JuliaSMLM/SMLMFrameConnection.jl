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
- `filter_high_chi2::Bool=false`: Deprecated and ignored (`true` logs a warning). It dropped
  tracks with a pair chi² above `chi2_filter_threshold` under the uncalibrated CRLB, which
  are the tracks that show k > 1, so it biased `k_scale` toward 1 (on SMLMSim data with
  true k 1.5 it gave 1.05). The fit itself leaves out mismatched pairs (see
  [`CalibrationResult`](@ref)).
- `chi2_filter_threshold::Float64=6.0`: Deprecated and ignored, with `filter_high_chi2`
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

The fit leaves out pairs whose chi² under the current fit (with the pair's full covariance,
σ_xy included) exceeds the chi²(2) 99.7% point (mismatched pairs), corrects the truncated bin means for the variance the cut removes, and
iterates to convergence, so `k_scale` and `sigma_motion_nm` are unbiased when the model holds.

# Fields
- `sigma_motion_nm::Float64`: Applied per-localization motion std dev in nm (√(A_applied/2) × 1000)
- `k_scale::Float64`: Applied CRLB scale factor (√B_applied)
- `A::Float64`: Intercept of the free two-parameter fit (pair-difference variance 2σ_motion², in μm²); can be negative
- `B::Float64`: Slope of the free two-parameter fit (k²). The line `A + B * bin_centers`
  goes through `bin_observed`. `B < 1` while `B_applied == 1` means the clamp bound
- `A_sigma::Float64`: Standard error of A
- `B_sigma::Float64`: Standard error of B
- `r_squared::Float64`: Weighted R² of the WLS fit (about the weighted mean)
- `mean_chi2::Float64`: Mean chi² (full pair covariance, uncalibrated CRLB) across all frame-to-frame pairs
- `n_pairs::Int`: Number of frame-to-frame pairs used in the final fit: kept by its chi² cut, inside the
  bins' span and in bins of 5 or more. On a fallback, the same count for the last fit
  attempted, or all pairs when there were fewer than 100
- `n_tracks_used::Int`: Number of tracks contributing pairs
- `n_tracks_filtered::Int`: Always 0 (the deprecated track filter is ignored)
- `bin_centers::Vector{Float64}`: Bin centers (CRLB variance) for diagnostic plots
- `bin_observed::Vector{Float64}`: Bin observed variance (truncation-corrected) of the free fit, for diagnostic plots
- `bin_counts::Vector{Int}`: Number of pairs per bin (parallel to bin_centers/bin_observed); use to grey/drop low-statistics bins in diagnostic plots
- `frame_shifts::Dict{Int, Vector{NTuple{2,Float64}}}`: Per-dataset (dx,dy) shifts for jitter plots
- `calibration_applied::Bool`: Whether calibration was actually applied (false on fallback)
- `warning::String`: Why calibration fell back (too few pairs or bins, or motion and CRLB
  scale not separable), or, on an applied calibration, that the fit did not converge
  (empty if OK)
- `A_applied::Float64`: Intercept applied, in μm² (≥ 0): the free fit's, or a refit with a
  bound held (A ≥ 0; B ≥ 1 with `clamp_k_to_one`). 0 on fallback
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

# Total order over every field of a localization, used to sort before the pipeline:
# localizations that tie on (dataset, frame, x, y) are ordered by their other fields,
# not by their place in the input; ones equal in every field are interchangeable.
_order_key(e) = (e.dataset, e.frame, e.x, e.y, e.σ_x, e.σ_y, e.σ_xy, e.photons, e.bg,
                 e.σ_photons, e.σ_bg, e.id, e.track_id)

# Covariance arithmetic. Every determinant in the package is built from these, so a
# covariance that passes the entry check (positive definite in exact arithmetic) keeps a
# positive determinant through every sum and inverse, however close |σ_xy| is to σ_xσ_y.

# Determinant σ_x²σ_y² - σ_xy² of a localization covariance, in Float64 as
# (p - q)(p + q) with p = σ_xσ_y carried exactly as a sum of two Float64 values (fma):
# the sign is exact for Float32 and Float64 input, and the factored form keeps the value
# accurate when q = |σ_xy| is close to p.
function _cov_det(σx, σy, σxy)
    a, b = Float64(σx), Float64(σy)
    p = a * b
    e = fma(a, b, -p)            # σ_xσ_y = p + e exactly
    q = abs(Float64(σxy))
    return ((p - q) + e) * ((p + q) + e)
end

# det(A + B) for 2×2 positive (semi)definite A = [a11 a12; a12 a22], B likewise, with
# their determinants da, db given accurately. det(A + B) = da + db + a11·b22 + a22·b11 -
# 2·a12·b12, and the cross term is written as a sum of non-negative parts,
# (√(a11·b22) - √(a22·b11))² + 2(√(a11·a22)·√(b11·b22) - a12·b12), the last one through
# √(a11·a22) - |a12| = da / (√(a11·a22) + |a12|): no subtraction of nearly equal values,
# so the result is positive whenever da or db is.
function _det_add(a11, a22, a12, da, b11, b22, b12, db)
    ra, rb = sqrt(a11 * a22), sqrt(b11 * b22)
    t1 = (sqrt(a11 * b22) - sqrt(a22 * b11))^2
    t2 = if a12 * b12 <= 0
        2 * (ra * rb + abs(a12 * b12))
    else
        2 * (da / (ra + abs(a12)) * rb + abs(a12) * db / (rb + abs(b12)))
    end
    return da + db + t1 + t2
end

# Determinant of the sum of two localization covariances (a pair difference's covariance).
function _pair_cov_det(σx1, σy1, σxy1, σx2, σy2, σxy2)
    return _det_add(Float64(σx1)^2, Float64(σy1)^2, Float64(σxy1), _cov_det(σx1, σy1, σxy1),
                    Float64(σx2)^2, Float64(σy2)^2, Float64(σxy2), _cov_det(σx2, σy2, σxy2))
end

# Precision-weighted mean of the localizations `idx` with full 2×2 covariances:
# Σ = (Σᵢ Σᵢ⁻¹)⁻¹, μ = Σ·Σᵢ Σᵢ⁻¹ xᵢ, in Float64. The determinant of the summed precision
# is carried along accurately (_det_add) instead of recomputed from the sums, where it
# cancels for nearly singular covariances. Each coordinate enters relative to an anchor:
# x relative to the localization with the largest x precision, y relative to the one
# with the largest y precision (first in `idx` on a tie). With a nearly singular
# precision, Σᵢ Σᵢ⁻¹ xᵢ cancels, and its error then scales with the offsets (about σ)
# instead of the absolute position, while each coordinate's dominant term stays exact.
# When every covariance in the group is the same, the result is exact in closed form,
# the arithmetic mean and Σ/n, however close to singular Σ is. The mean is taken as
# x₁ + Σ (xᵢ - x₁)/n: exactly x₁ when all coordinates are equal (a mean that rounds
# away from them reads as a real separation downstream, in the density estimate; x₁
# itself, so a group of one, or of equal coordinates, keeps even the sign of a zero), and
# each difference is scaled before the sum, so it overflows only when two coordinates
# differ by more than the largest Float64.
# Returns (μx, μy, Σ11, Σ22, Σ12).
function _precision_mean(x, y, σx, σy, σxy, idx)
    i1 = first(idx)
    if all(i -> σx[i] == σx[i1] && σy[i] == σy[i1] && σxy[i] == σxy[i1], idx)
        n = length(idx)
        x1, y1 = Float64(x[i1]), Float64(y[i1])
        sx = sum((Float64(x[i]) - x1) / n for i in idx)
        sy = sum((Float64(y[i]) - y1) / n for i in idx)
        return sx == 0 ? x1 : x1 + sx, sy == 0 ? y1 : y1 + sy,
               Float64(σx[i1])^2 / n, Float64(σy[i1])^2 / n, Float64(σxy[i1]) / n
    end
    ix, iy, wx, wy = i1, i1, -Inf, -Inf
    for i in idx
        d = _cov_det(σx[i], σy[i], σxy[i])
        px, py = Float64(σy[i])^2 / d, Float64(σx[i])^2 / d   # P11, P22
        px > wx && ((ix, wx) = (i, px))
        py > wy && ((iy, wy) = (i, py))
    end
    S11 = S22 = S12 = det_S = 0.0
    m1 = m2 = 0.0
    x0, y0 = Float64(x[ix]), Float64(y[iy])
    for (j, i) in enumerate(idx)
        det_i = _cov_det(σx[i], σy[i], σxy[i])   # > 0: the entry check holds
        # Precision matrix P = [σ_y² -σ_xy; -σ_xy σ_x²] / det, with det(P) = 1/det
        P11 = Float64(σy[i])^2 / det_i
        P22 = Float64(σx[i])^2 / det_i
        P12 = -Float64(σxy[i]) / det_i
        det_S = j == 1 ? 1 / det_i : _det_add(S11, S22, S12, det_S, P11, P22, P12, 1 / det_i)
        S11 += P11
        S22 += P22
        S12 += P12
        dx, dy = Float64(x[i]) - x0, Float64(y[i]) - y0
        m1 += P11 * dx + P12 * dy
        m2 += P12 * dx + P22 * dy
    end
    Σ11, Σ22, Σ12 = S22 / det_S, S11 / det_S, -S12 / det_S
    return x0 + (Σ11 * m1 + Σ12 * m2), y0 + (Σ12 * m1 + Σ22 * m2), Σ11, Σ22, Σ12
end

# Round a positive-definite covariance (σx, σy, σxy, computed in Float64) to the element
# type ET. Rounding can carry a correlation within a rounding error of ±1 past it; then
# σxy becomes the largest |σxy| in ET below σx·σy, a few steps from that bound at most
# (and 0 if those steps do not reach it, so this always returns).
function _round_cov(::Type{ET}, σx, σy, σxy) where {ET}
    ox, oy, oxy = ET(σx), ET(σy), ET(σxy)
    _cov_det(ox, oy, oxy) > 0 && return ox, oy, oxy
    (ox > 0 && oy > 0 && isfinite(ox) && isfinite(oy)) || return ox, oy, oxy
    c = ET(copysign(Float64(ox) * Float64(oy), Float64(σxy)))
    for _ in 1:8
        _cov_det(ox, oy, c) > 0 && return ox, oy, c
        c = c > 0 ? prevfloat(c) : nextfloat(c)
    end
    return ox, oy, zero(ET)
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
        (isfinite(e.σ_xy) && _cov_det(e.σ_x, e.σ_y, e.σ_xy) > 0) ||
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
