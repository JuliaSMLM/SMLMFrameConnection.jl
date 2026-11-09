# Uncertainty calibration for frame-connected tracks.
#
# Analyzes frame-to-frame jitter within connected tracks to estimate motion
# variance (σ_motion²) and CRLB scale factor (k²), then applies corrected
# uncertainties: Σ_corrected = σ_motion² I + k² Σ_CRLB

# Pairs whose chi² under the current fit exceeds the chi²(2) 99.7% point are left out
# of the fit as mismatches.
const _CAL_CHI2_CUT = -2 * log(0.003)
# The cut statistic is S = dᵀ S_pair⁻¹ d with the full pair covariance S_pair (x–y term
# included). Whitened, d = L·w with w ~ N(0, I) and S = |w|², so
# E[d dᵀ | S ≤ c] = L·E[w wᵀ | S ≤ c]·Lᵀ = g(c)·S_pair, with
# g(c) = E[S | S ≤ c] / 2 = 1 - (c/2)·e^(-c/2) / (1 - e^(-c/2)) for S ~ chi²(2).
# Dividing the truncated per-axis mean square by g(c) makes it an unbiased estimate of
# the per-axis pair variance again, whatever the correlation.
const _CAL_TRUNC_FACTOR = 1 - (_CAL_CHI2_CUT / 2) * exp(-_CAL_CHI2_CUT / 2) /
                              (1 - exp(-_CAL_CHI2_CUT / 2))
const _CAL_MAX_ITER = 100
# The fitted slope B = k² is used only when B ≥ this many standard errors.
const _CAL_MIN_SLOPE_SIGNIF = 3.0
# Fewest frame-to-frame pairs a calibration is applied from (after the cut).
const _CAL_MIN_PAIRS = 100
# Reported precisions whose range is at most this fraction of the largest are one value
# (one bin): equal to within rounding.
const _CAL_EQUAL_REL = 1e-12

# chi² of a pair difference (dx, dy) under the pair covariance [sxx sxy; sxy syy] with
# determinant `det` (given, computed without cancellation); it is chi²(2) when that
# covariance is right. Written as dx²/sxx + (dy - (sxy/sxx)·dx)²·sxx/det, a sum of
# squares. Inf when the covariance is not positive definite: the cut then leaves the
# pair out.
function _pair_chi2(dx, dy, sxx, syy, sxy, det)
    (sxx > 0 && syy > 0 && det > 0) || return Inf
    return dx^2 / sxx + sxx * (dy - (sxy / sxx) * dx)^2 / det
end

# Pairs kept by the cut: chi² under the model covariance A·I + B·(Σ₁ + Σ₂), whose
# determinant A² + A·B·(vx + vy) + B²·det(Σ₁ + Σ₂) is a sum of non-negative terms (A ≥ 0).
function _cut_keep(dx, dy, vx, vy, cxy, dS, A::Float64, B::Float64)
    keep = falses(length(dx))
    @inbounds for i in eachindex(dx)
        detC = A^2 + A * B * (vx[i] + vy[i]) + B^2 * dS[i]
        keep[i] = _pair_chi2(dx[i], dy[i], A + B * vx[i], A + B * vy[i], B * cxy[i], detC) <=
                  _CAL_CHI2_CUT
    end
    return keep
end

# Weighted least squares of binned `observed` on binned `reported` over the pairs in
# `keep`; `bin_of[i]` is pair i's bin. Bins with fewer than 5 pairs are skipped. Each bin
# mean is weighted by n / v², with v = A + B * center the variance the current fit (A, B)
# predicts: a pair's (dx² + dy²)/2 is v·chi²(2)/2, whose variance is v². Weighting by a
# bin's own sample variance instead favors bins whose mean came out low (for chi²-like
# values the sample mean and variance rise together), which biased k low on small data.
# `fixA` or `fixB` holds that parameter at the given value and fits the other alone (the
# refits of a determined free fit whose bound binds).
# Every fit needs 3 bins: the two-parameter fit one beyond its parameters, for a standard
# error, and a bound refit the same support as the free fit it replaces (its own cut can
# empty bins; a refit left with fewer has not been checked against the range).
# Returns a NamedTuple (with `n_used`, the pairs in the bins the fit used), or a String
# naming why no fit was possible.
function _binned_wls(reported, observed, bin_of, keep, n_bins, A, B;
                     fixA::Union{Nothing,Float64} = nothing,
                     fixB::Union{Nothing,Float64} = nothing)
    bin_centers = Float64[]
    bin_observed = Float64[]
    bin_var = Float64[]
    n_per_bin = Int[]
    members = [Int[] for _ in 1:n_bins]
    for i in eachindex(reported)
        keep[i] && push!(members[bin_of[i]], i)
    end
    for m in members
        if length(m) >= 5
            push!(bin_centers, mean(@view reported[m]))
            push!(bin_observed, mean(@view observed[m]))
            push!(bin_var, var(@view observed[m]))
            push!(n_per_bin, length(m))
        end
    end

    n_valid_bins = length(bin_centers)
    n_valid_bins < 3 && return "Too few valid bins ($n_valid_bins < 3)"

    # WLS on binned data: weight = 1/var_of_bin_mean (model variance), normalized
    # A negative intercept can make the prediction small or negative in the lowest bins;
    # floor it at half the slope term there (still model-based, not the bin's own
    # sample variance). Only a non-positive slope falls back to the sample variance.
    var_of_bin_mean = if B > 0
        predicted = max.(A .+ B .* bin_centers, 0.5 .* B .* bin_centers)
        predicted .^ 2 ./ n_per_bin
    else
        bin_var ./ n_per_bin
    end
    w_model = [v > 0 ? 1.0 / v : 0.0 for v in var_of_bin_mean]
    sum(w_model) > 0 || return "All bin weights are zero"
    weights = w_model ./ sum(w_model) .* n_valid_bins

    # Weighted R²: residual and total sums of squares both about the weighted fit, the
    # total about the weighted mean.
    weighted_mean = sum(weights .* bin_observed) / sum(weights)
    ss_tot = sum(weights .* (bin_observed .- weighted_mean) .^ 2)
    if fixA === nothing && fixB === nothing
        # Centered at the weighted mean of the bin centers. The normal equations of
        # (1, x) cancel when the reported precisions barely vary: their determinant is
        # Σw·Sxx formed as a difference of two nearly equal products, which rounded to
        # a nonzero value of a singular matrix and made the solve throw.
        sw = sum(weights)
        xbar = sum(weights .* bin_centers) / sw
        ybar = sum(weights .* bin_observed) / sw
        dxc = bin_centers .- xbar
        Sxx = sum(weights .* dxc .^ 2)
        # Precisions equal to within rounding share one bin (_calibration_bins), so
        # distinct bins differ; this guards only the impossible exact zero.
        Sxx > 0 || return "Reported precisions do not vary across bins"
        Bfit = sum(weights .* dxc .* (bin_observed .- ybar)) / Sxx
        Afit = ybar - Bfit * xbar
        ss_res = sum(weights .* (bin_observed .- (Afit .+ Bfit .* bin_centers)) .^ 2)
        s2 = ss_res / max(1, n_valid_bins - 2)
        A_sigma = sqrt(max(0.0, s2 * (1 / sw + xbar^2 / Sxx)))
        B_sigma = sqrt(max(0.0, s2 / Sxx))
    elseif fixA !== nothing
        Afit = fixA
        swcc = sum(weights .* bin_centers .^ 2)
        Bfit = sum(weights .* bin_centers .* (bin_observed .- Afit)) / swcc
        ss_res = sum(weights .* (bin_observed .- (Afit .+ Bfit .* bin_centers)) .^ 2)
        A_sigma = 0.0
        B_sigma = sqrt(max(0.0, ss_res / (n_valid_bins - 1) / swcc))
    else
        Bfit = fixB
        Afit = sum(weights .* (bin_observed .- Bfit .* bin_centers)) / sum(weights)
        ss_res = sum(weights .* (bin_observed .- (Afit .+ Bfit .* bin_centers)) .^ 2)
        A_sigma = sqrt(max(0.0, ss_res / (n_valid_bins - 1) / sum(weights)))
        B_sigma = 0.0
    end
    r_squared = ss_tot > 0 ? 1.0 - ss_res / ss_tot : NaN

    return (A = Afit, B = Bfit, A_sigma = A_sigma, B_sigma = B_sigma,
            r_squared = r_squared, bin_centers = bin_centers,
            bin_observed = bin_observed, n_per_bin = n_per_bin, n_used = sum(n_per_bin))
end

# Equal-width bins over all pairs' reported values, fixed across iterations: bin b holds
# x_min + (b-1)·w ≤ r < x_min + b·w, decided by comparison with those edges (not by
# rounding (r - x_min)/w, which can put a value on an edge in the bin below), and the last
# bin everything above its lower edge. Values whose range is at most _CAL_EQUAL_REL of the
# largest are one value in one bin: split one ulp apart, their bins carry no slope.
function _calibration_bins(reported, n_bins)
    x_min, x_max = extrema(reported)
    x_max - x_min > _CAL_EQUAL_REL * x_max || return ones(Int, length(reported))
    w = (x_max - x_min) / n_bins
    edges = [x_min + b * w for b in 1:n_bins-1]
    return [searchsortedlast(edges, r) + 1 for r in reported]
end

# Iterate cut, truncation-corrected binned fit, and weights from (A, B) to a fixed
# point; `fixA`/`fixB` as in `_binned_wls`. `dx`, `dy` are the pair differences, `vx`,
# `vy`, `cxy` the sums of the pair's reported variances and covariances, `dS` the
# determinant of that sum (_pair_cov_det). Returns (fit, keep, converged).
function _iterate_calibration(reported, observed, bin_of, n_bins, dx, dy, vx, vy, cxy, dS,
                              A, B; fixA = nothing, fixB = nothing)
    fit = "No calibration fit was attempted"
    keep = trues(length(reported))
    converged = false
    for _ in 1:_CAL_MAX_ITER
        # The cut uses the pair's full model covariance A·I + B·(Σ₁ + Σ₂), x–y term
        # included: only then is its statistic chi²(2) and the truncation factor right.
        # A trial A < 0 (σ_motion² < 0, outside the model) enters the cut as 0: with a
        # strongly correlated Σ, A·I + B·Σ is not positive definite for any A below
        # -B·λ_min(Σ), and every such pair would be cut. Its determinant,
        # A² + A·B·(vx + vy) + B²·det(Σ₁ + Σ₂), is a sum of non-negative terms.
        Ac = max(A, 0.0)
        keep = _cut_keep(dx, dy, vx, vy, cxy, dS, Ac, B)
        fit = _binned_wls(reported, observed, bin_of, keep, n_bins, A, B;
                          fixA = fixA, fixB = fixB)
        fit isa String && break
        done = abs(fit.B - B) <= 1e-6 * max(1.0, abs(B)) &&
               abs(fit.A - A) <= 1e-6 * max(1.0, mean(reported))
        A, B = fit.A, fit.B
        converged = done
        done && break
    end
    return fit, keep, converged
end

"""
    analyze_calibration(smld::BasicSMLD, config::CalibrationConfig) -> CalibrationResult

Analyze frame-to-frame jitter in connected tracks to estimate uncertainty
calibration parameters.

Collects consecutive-frame pairs within each track, computes observed positional
variance vs CRLB variance, and fits the model: `observed_var = A + B * CRLB_var`
using weighted least squares on binned data. `A = 2σ_motion²` (both localizations of a
pair carry motion jitter) and `B = k²`. Pairs whose chi² under the current fit (with the
pair's full covariance, σ_xy included) exceeds the chi²(2) 99.7% point are left out as
mismatches; the truncated bin means are corrected for the variance the cut removes, and
the fit iterates to convergence.

The fit is applied only when it separates the two terms: slope positive and at least 3
standard errors, and R² ≥ 0.1. A determined fit whose bound binds (A < 0, or B < 1 with
`clamp_k_to_one`) is refit with that parameter held at its bound.

Falls back gracefully (returns `calibration_applied=false`, with the reason in `warning`)
if fewer than 100 frame-to-frame pairs are available, the reported precisions fill fewer
than 3 bins (reported precisions whose range is at most 1e-12 of the largest are one
value, one bin), the fit cannot separate the two terms (a narrow range of reported
precisions, or a dominant motion term), a bound refit keeps fewer than 3 bins after its
own cut, or fewer than 100 pairs remain in the applied fit after the cut.
"""
function analyze_calibration(smld::BasicSMLD{T,E},
                             config::CalibrationConfig) where {T, E<:SMLMData.AbstractEmitter}
    emitters = smld.emitters
    n = length(emitters)

    # Group emitters by track_id
    track_map = Dict{Int, Vector{Int}}()
    for i in 1:n
        tid = emitters[i].track_id
        tid == 0 && continue
        if haskey(track_map, tid)
            push!(track_map[tid], i)
        else
            track_map[tid] = [i]
        end
    end

    # Collect frame-to-frame pairs and frame shifts
    dx_vec = Float64[]
    dy_vec = Float64[]
    var_x_vec = Float64[]  # sum of CRLB variances for pair
    var_y_vec = Float64[]
    cov_xy_vec = Float64[]  # sum of CRLB covariances σ_xy for pair
    det_vec = Float64[]     # determinant of the pair's summed CRLB covariance
    n_tracks_used = 0
    frame_shifts = Dict{Int, Vector{NTuple{2,Float64}}}()

    for (tid, indices) in track_map
        length(indices) < 2 && continue

        # Sort by frame within track
        sorted = sort(indices, by=i -> emitters[i].frame)
        n_before = length(dx_vec)

        for j in 1:(length(sorted) - 1)
            i1, i2 = sorted[j], sorted[j+1]
            e1, e2 = emitters[i1], emitters[i2]

            # Only consecutive frames in same dataset
            e1.dataset != e2.dataset && continue
            abs(e2.frame - e1.frame) != 1 && continue

            # Each field promoted to Float64 before any arithmetic: squares and sums
            # formed in Float32 can turn two positive-definite covariances into a pair
            # covariance that is not.
            dx = Float64(e2.x) - Float64(e1.x)
            dy = Float64(e2.y) - Float64(e1.y)
            vx = Float64(e1.σ_x)^2 + Float64(e2.σ_x)^2
            vy = Float64(e1.σ_y)^2 + Float64(e2.σ_y)^2
            cxy = Float64(e1.σ_xy) + Float64(e2.σ_xy)

            push!(dx_vec, dx)
            push!(dy_vec, dy)
            push!(var_x_vec, vx)
            push!(var_y_vec, vy)
            push!(cov_xy_vec, cxy)
            push!(det_vec, _pair_cov_det(e1.σ_x, e1.σ_y, e1.σ_xy, e2.σ_x, e2.σ_y, e2.σ_xy))

            # Store frame shifts for jitter diagnostic
            ds = e1.dataset
            shift = (dx, dy)
            if haskey(frame_shifts, ds)
                push!(frame_shifts[ds], shift)
            else
                frame_shifts[ds] = [shift]
            end
        end

        length(dx_vec) > n_before && (n_tracks_used += 1)
    end

    # The track filter is deprecated and ignored (see CalibrationConfig): it dropped tracks
    # by their chi² under the uncalibrated CRLB, which are the tracks that show k > 1, and
    # biased k_scale toward 1 (true k 1.5 gave 1.05). The fit's own cut leaves out
    # mismatched pairs. Warned before any return, so small data hear it too.
    config.filter_high_chi2 && @warn "CalibrationConfig(filter_high_chi2 = true) is " *
        "deprecated and ignored: dropping tracks by their uncalibrated chi² biased " *
        "k_scale toward 1. The fit already leaves out mismatched pairs." maxlog = 1

    n_pairs = length(dx_vec)
    mean_chi2 = n_pairs > 0 ?
        mean(_pair_chi2.(dx_vec, dy_vec, var_x_vec, var_y_vec, cov_xy_vec, det_vec)) : NaN

    # Fallback: not enough data
    if n_pairs < _CAL_MIN_PAIRS
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Too few frame-to-frame pairs ($n_pairs < $_CAL_MIN_PAIRS)")
    end

    # Model, per axis, for the difference of two consecutive localizations of one
    # emitter: Var(Δx) = A + B·(σ_x1² + σ_x2²). Each localization carries its own
    # motion jitter, so the intercept is A = 2·σ_motion²; the slope is B = k².
    #
    # Fit robustly to mismatched pairs (double fits, wrong links) without biasing the
    # scale: a pair is left out when its chi² under the CURRENT fit exceeds the
    # chi²(2) 99.7% point, the truncated bin means are divided by the factor the
    # truncation removes, and the fit repeats until it settles. A fixed cut on the
    # uncalibrated chi² (the earlier method) removes exactly the large pairs that
    # carry the excess when the CRLB is too small, and biases k toward 1.
    vx_nm2 = var_x_vec .* 1e6
    vy_nm2 = var_y_vec .* 1e6
    cxy_nm2 = cov_xy_vec .* 1e6
    det_nm4 = det_vec .* 1e12
    dx_nm = dx_vec .* 1e3
    dy_nm = dy_vec .* 1e3
    # Reported value per pair and its bin, in the same arithmetic as before this
    # calibration was reworked: equal values must land in one bin, and a value on a bin
    # edge in the same bin as before (another order of operations can round two equal
    # sums apart, or across an edge, and drop a bin).
    reported_all = ((var_x_vec .+ var_y_vec) ./ 2) .* 1e6
    observed_all = (((dx_vec .^ 2 .+ dy_vec .^ 2) ./ 2) .* 1e6) ./ _CAL_TRUNC_FACTOR

    n_bins = 20
    bin_of = _calibration_bins(reported_all, n_bins)

    # Start from the median: with A = 0, chi²/B is chi²(2), whose median is 2·log(2).
    z_raw = _pair_chi2.(dx_nm, dy_nm, vx_nm2, vy_nm2, cxy_nm2, det_nm4)
    A = 0.0
    B = median(z_raw) / (2 * log(2))

    iter(A0, B0; kw...) = _iterate_calibration(reported_all, observed_all, bin_of, n_bins,
                                               dx_nm, dy_nm, vx_nm2, vy_nm2, cxy_nm2,
                                               det_nm4, A0, B0; kw...)
    fallback(msg, n) = _fallback_result(mean_chi2, n, n_tracks_used, frame_shifts, msg)
    too_few(n) = "Too few pairs in the fit after the chi² cut ($n < $_CAL_MIN_PAIRS)"
    # Pairs a fit actually used: kept by the cut and in a bin with enough pairs.
    used(f) = f isa String ? 0 : f.n_used

    fit, keep, converged = iter(A, B)
    # No fit when the reported precisions fill fewer than 3 bins (a narrow or single
    # value).
    fit isa String && return fallback(fit, count(keep))
    free = fit   # reported as A, B, errors, R² and the bins

    # The free fit separates a motion term from the CRLB scale only when its slope is
    # positive and determined (B > 0, B ≥ 3 standard errors, R² ≥ 0.1). When it is not,
    # the reported precisions span too narrow a range to tell the two terms apart, or a
    # motion term dominates and the observed variance barely changes with the CRLB: the
    # calibration falls back rather than apply an undetermined split.
    determined = free.B > 0 && free.B >= _CAL_MIN_SLOPE_SIGNIF * free.B_sigma &&
                 free.r_squared >= 0.1
    determined || return fallback("Motion and CRLB scale not separable (free fit: B = " *
        "$(round(free.B, sigdigits=3)) ± $(round(free.B_sigma, sigdigits=3)), " *
        "R² = $(round(free.r_squared, digits=3)))", used(free))
    # A determined free fit with A < 0 refits the scale alone: the bound σ_motion² ≥ 0
    # binds, and keeping the slope fitted jointly with the negative intercept moves the
    # applied variance away from the data.
    if fit.A < 0
        fit, keep, conv = iter(0.0, fit.B; fixA = 0.0)
        converged &= conv
        fit isa String && return fallback("Refit with σ_motion = 0: $fit", count(keep))
    end
    if config.clamp_k_to_one && fit.B < 1
        fit, keep, conv = iter(max(fit.A, 0.0), 1.0; fixB = 1.0)
        converged &= conv
        fit isa String && return fallback("Refit with k = 1: $fit", count(keep))
        if fit.A < 0   # both bounds bind: the CRLB as reported, no motion term
            fit = merge(fit, (A = 0.0, A_sigma = 0.0))
        end
    end
    # Every refit runs its own cut: the applied fit, not only the free one, must use
    # enough pairs (kept by its cut and in its bins).
    n_pairs = used(fit)
    n_pairs < _CAL_MIN_PAIRS && return fallback(too_few(n_pairs), n_pairs)
    # Applied: within the bounds (A ≥ 0; B ≥ 1 with the clamp). Fit was in nm².
    A_applied, B_applied = fit.A, fit.B
    sigma_motion_nm = sqrt(A_applied / 2)   # A is the pair intercept 2·σ_motion²
    k_scale = sqrt(B_applied)                # B is k²

    # Reported fit: the free two-parameter fit (A, B, errors, R², bins), so the line
    # A + B·x goes through the stored bins. Converted from nm² to μm².
    (; A, B, A_sigma, B_sigma, r_squared, bin_centers, bin_observed, n_per_bin) = free

    return CalibrationResult(
        sigma_motion_nm, k_scale,
        A * 1e-6, B, A_sigma * 1e-6, B_sigma, r_squared, mean_chi2,
        n_pairs, n_tracks_used, 0,   # n_tracks_filtered: the track filter is ignored
        bin_centers .* 1e-6, bin_observed .* 1e-6, copy(n_per_bin), frame_shifts,
        true, converged ? "" :
              "Calibration fit did not converge in $(_CAL_MAX_ITER) iterations; the last fit was applied",
        A_applied * 1e-6, B_applied
    )
end

"""
    apply_calibration(smld::BasicSMLD, result::CalibrationResult) -> BasicSMLD

Apply calibrated uncertainties to emitters. Returns a new SMLD with corrected
covariance: `Σ_corrected = σ_motion² I + k² Σ_CRLB`.

Only modifies σ_x, σ_y, σ_xy. Positions and other fields are unchanged.
"""
function apply_calibration(smld::BasicSMLD{T,E},
                           result::CalibrationResult) where {T, E<:SMLMData.AbstractEmitter}
    !result.calibration_applied && return smld

    sigma_motion_sq = (result.sigma_motion_nm / 1000)^2  # nm to μm, then squared
    k_sq = result.k_scale^2

    emitters = smld.emitters
    ET = typeof(first(emitters).x)
    new_emitters = Vector{SMLMData.Emitter2DFit{ET}}(undef, length(emitters))

    for i in eachindex(emitters)
        e = emitters[i]
        # Σ_corrected = σ_motion² I + k² Σ_CRLB
        # In Float64 (squares of tiny Float32 values underflow), then rounded to ET. A
        # positive-definite input stays positive definite: _round_cov keeps it so through
        # the rounding. Other input (only reachable by calling this directly;
        # frameconnect refuses it) is converted as computed.
        sx = sqrt(sigma_motion_sq + k_sq * Float64(e.σ_x)^2)
        sy = sqrt(sigma_motion_sq + k_sq * Float64(e.σ_y)^2)
        sxy = k_sq * Float64(e.σ_xy)
        new_σx, new_σy, new_σxy = _cov_det(e.σ_x, e.σ_y, e.σ_xy) > 0 ?
            _round_cov(ET, sx, sy, sxy) : (ET(sx), ET(sy), ET(sxy))

        new_emitters[i] = SMLMData.Emitter2DFit{ET}(
            e.x, e.y, e.photons, e.bg,
            new_σx, new_σy, new_σxy,
            e.σ_photons, e.σ_bg,
            e.frame, e.dataset, e.track_id, e.id
        )
    end

    return BasicSMLD(new_emitters, smld.camera, smld.n_frames,
                     smld.n_datasets, copy(smld.metadata))
end

function _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts, warning)
    return CalibrationResult(
        0.0, 1.0,               # sigma_motion_nm, k_scale
        0.0, 1.0, 0.0, 0.0,    # A, B, A_sigma, B_sigma
        0.0, mean_chi2,         # r_squared, mean_chi2
        n_pairs, n_tracks_used, 0,   # n_tracks_filtered: the track filter is ignored
        Float64[], Float64[], Int[],   # bin_centers, bin_observed, bin_counts
        frame_shifts,
        false, warning,
        0.0, 1.0                # A_applied, B_applied: no correction
    )
end
