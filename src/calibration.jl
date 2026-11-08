# Uncertainty calibration for frame-connected tracks.
#
# Analyzes frame-to-frame jitter within connected tracks to estimate motion
# variance (σ_motion²) and CRLB scale factor (k²), then applies corrected
# uncertainties: Σ_corrected = σ_motion² I + k² Σ_CRLB

# Pairs whose chi² under the current fit exceeds the chi²(2) 99.7% point are left out
# of the fit as mismatches.
const _CAL_CHI2_CUT = -2 * log(0.003)
# For S ~ chi²(2) = u₁ + u₂ (u iid chi²(1)), E[u₁ | S ≤ c] = E[S | S ≤ c] / 2
# = 1 - (c/2)·e^(-c/2) / (1 - e^(-c/2)). Dividing the truncated per-axis mean square by
# this factor makes it an unbiased estimate of the per-axis pair variance again.
const _CAL_TRUNC_FACTOR = 1 - (_CAL_CHI2_CUT / 2) * exp(-_CAL_CHI2_CUT / 2) /
                              (1 - exp(-_CAL_CHI2_CUT / 2))
const _CAL_MAX_ITER = 100
# The fitted slope B = k² is used only when B ≥ this many standard errors.
const _CAL_MIN_SLOPE_SIGNIF = 3.0

# Weighted least squares of binned `observed` on binned `reported` over the pairs in
# `keep`; `bin_of[i]` is pair i's bin. Bins with fewer than 5 pairs are skipped. Each bin
# mean is weighted by n / v², with v = A + B * center the variance the current fit (A, B)
# predicts: a pair's (dx² + dy²)/2 is v·chi²(2)/2, whose variance is v². Weighting by a
# bin's own sample variance instead favors bins whose mean came out low (for chi²-like
# values the sample mean and variance rise together), which biased k low on small data.
# `fixA` or `fixB` holds that parameter at the given value and fits the other alone.
# Returns a NamedTuple, or a String naming why no fit was possible.
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
    weights = [v > 0 ? 1.0 / v : 0.0 for v in var_of_bin_mean]
    sum(weights) ≈ 0 && return "All bin weights are zero"
    weights = weights ./ sum(weights) .* n_valid_bins

    ss_tot = sum(weights .* (bin_observed .- mean(bin_observed)) .^ 2)
    if fixA === nothing && fixB === nothing
        X = hcat(ones(n_valid_bins), bin_centers)
        W = Diagonal(weights)
        XtWX = X' * W * X
        XtWy = X' * W * bin_observed
        det_XtWX = XtWX[1,1] * XtWX[2,2] - XtWX[1,2]^2
        abs(det_XtWX) < eps() && return "Singular WLS matrix"
        coeffs = XtWX \ XtWy
        Afit, Bfit = coeffs[1], coeffs[2]
        ss_res = sum(weights .* (bin_observed .- (Afit .+ Bfit .* bin_centers)) .^ 2)
        var_coeffs = ss_res / max(1, n_valid_bins - 2) * inv(XtWX)
        A_sigma = sqrt(max(0.0, var_coeffs[1, 1]))
        B_sigma = sqrt(max(0.0, var_coeffs[2, 2]))
    elseif fixA !== nothing
        Afit = fixA
        swcc = sum(weights .* bin_centers .^ 2)
        Bfit = sum(weights .* bin_centers .* (bin_observed .- Afit)) / swcc
        ss_res = sum(weights .* (bin_observed .- (Afit .+ Bfit .* bin_centers)) .^ 2)
        A_sigma = 0.0
        B_sigma = sqrt(max(0.0, ss_res / max(1, n_valid_bins - 1) / swcc))
    else
        Bfit = fixB
        Afit = sum(weights .* (bin_observed .- Bfit .* bin_centers)) / sum(weights)
        ss_res = sum(weights .* (bin_observed .- (Afit .+ Bfit .* bin_centers)) .^ 2)
        A_sigma = sqrt(max(0.0, ss_res / max(1, n_valid_bins - 1) / sum(weights)))
        B_sigma = 0.0
    end
    r_squared = ss_tot > 0 ? 1.0 - ss_res / ss_tot : NaN
    # Misfit comparable between models: sum over bins of n·(observed/predicted - 1)²
    # (each bin's mean has relative variance about 1/n under the model), per degree of
    # freedom left after the fitted parameters.
    npar = (fixA === nothing) + (fixB === nothing)
    pred = max.(Afit .+ Bfit .* bin_centers, eps())
    misfit = sum(n_per_bin .* (bin_observed ./ pred .- 1) .^ 2) / max(1, n_valid_bins - npar)

    return (A = Afit, B = Bfit, A_sigma = A_sigma, B_sigma = B_sigma,
            r_squared = r_squared, misfit = misfit, bin_centers = bin_centers,
            bin_observed = bin_observed, n_per_bin = n_per_bin)
end

# Iterate cut, truncation-corrected binned fit, and weights from (A, B) to a fixed
# point; `fixA`/`fixB` as in `_binned_wls`. Returns (fit, keep, converged).
function _iterate_calibration(reported, observed, bin_of, n_bins, vx, vy, dx2, dy2, A, B;
                              fixA = nothing, fixB = nothing)
    fit = "No calibration fit was attempted"
    keep = trues(length(reported))
    converged = false
    for _ in 1:_CAL_MAX_ITER
        sx = max.(A .+ B .* vx, eps())
        sy = max.(A .+ B .* vy, eps())
        keep = (dx2 ./ sx .+ dy2 ./ sy) .<= _CAL_CHI2_CUT
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
pair carry motion jitter) and `B = k²`. Pairs whose chi² under the current fit exceeds
the chi²(2) 99.7% point are left out as mismatches; the truncated bin means are
corrected for the variance the cut removes, and the fit iterates to convergence.

When the free fit cannot separate the two terms (slope not positive and at least 3
standard errors, or R² < 0.1: a narrow range of reported precisions, or a dominant
motion term), both one-parameter models are fit, the CRLB scale alone (σ_motion = 0)
and the motion term alone (k = 1), the better-fitting one is applied, and `warning`
says which (also as a logged warning). When both a motion term and a CRLB scale matter
but the data cannot separate them (typically 800 pairs or fewer), neither
one-parameter model is exact: the applied variance stays close near the median σ but
can be off by tens of percent at the extremes. Falls back gracefully (returns
`calibration_applied=false`) if fewer than 100 frame-to-frame pairs are available or
neither one-parameter model is determined.
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
    pair_tid_vec = Int[]    # track ID for each pair (for filtering)
    chi2_per_track = Dict{Int, Vector{Float64}}()
    frame_shifts = Dict{Int, Vector{NTuple{2,Float64}}}()

    for (tid, indices) in track_map
        length(indices) < 2 && continue

        # Sort by frame within track
        sorted = sort(indices, by=i -> emitters[i].frame)
        track_chi2s = Float64[]

        for j in 1:(length(sorted) - 1)
            i1, i2 = sorted[j], sorted[j+1]
            e1, e2 = emitters[i1], emitters[i2]

            # Only consecutive frames in same dataset
            e1.dataset != e2.dataset && continue
            abs(e2.frame - e1.frame) != 1 && continue

            dx = Float64(e2.x - e1.x)
            dy = Float64(e2.y - e1.y)
            vx = Float64(e1.σ_x^2 + e2.σ_x^2)
            vy = Float64(e1.σ_y^2 + e2.σ_y^2)

            chi2 = dx^2 / max(vx, eps()) + dy^2 / max(vy, eps())

            push!(dx_vec, dx)
            push!(dy_vec, dy)
            push!(var_x_vec, vx)
            push!(var_y_vec, vy)
            push!(pair_tid_vec, tid)
            push!(track_chi2s, chi2)

            # Store frame shifts for jitter diagnostic
            ds = e1.dataset
            shift = (dx, dy)
            if haskey(frame_shifts, ds)
                push!(frame_shifts[ds], shift)
            else
                frame_shifts[ds] = [shift]
            end
        end

        if !isempty(track_chi2s)
            chi2_per_track[tid] = track_chi2s
        end
    end

    n_pairs = length(dx_vec)
    n_tracks_used = length(chi2_per_track)
    mean_chi2 = n_pairs > 0 ? mean(dx_vec .^ 2 ./ max.(var_x_vec, eps()) .+
                                    dy_vec .^ 2 ./ max.(var_y_vec, eps())) : NaN

    # Fallback: not enough data
    if n_pairs < 100
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Too few frame-to-frame pairs ($n_pairs < 100)")
    end

    # Optional chi² track filtering
    n_tracks_filtered = 0
    if config.filter_high_chi2
        # Identify tracks with any chi² above threshold
        filtered_tids = Set{Int}()
        for (tid, chi2s) in chi2_per_track
            if any(c -> c > config.chi2_filter_threshold, chi2s)
                push!(filtered_tids, tid)
            end
        end
        n_tracks_filtered = length(filtered_tids)

        # Build keep mask using stored per-pair track IDs
        keep_mask = [!(tid in filtered_tids) for tid in pair_tid_vec]

        dx_vec = dx_vec[keep_mask]
        dy_vec = dy_vec[keep_mask]
        var_x_vec = var_x_vec[keep_mask]
        var_y_vec = var_y_vec[keep_mask]
        n_pairs = length(dx_vec)

        if n_pairs < 100
            return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                    "Too few pairs after chi² filtering ($n_pairs < 100)";
                                    n_tracks_filtered=n_tracks_filtered)
        end
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
    dx2_nm2 = dx_vec .^ 2 .* 1e6
    dy2_nm2 = dy_vec .^ 2 .* 1e6
    reported_all = (vx_nm2 .+ vy_nm2) ./ 2
    observed_all = ((dx2_nm2 .+ dy2_nm2) ./ 2) ./ _CAL_TRUNC_FACTOR

    # Equal-width bin edges over all pairs, fixed across iterations.
    x_min, x_max = extrema(reported_all)
    n_bins = 20
    bin_width = (x_max - x_min) / n_bins
    bin_of = [bin_width > 0 ? clamp(floor(Int, (r - x_min) / bin_width) + 1, 1, n_bins) : 1
              for r in reported_all]

    # Start from the median: with A = 0, chi²/B is chi²(2), whose median is 2·log(2).
    z_raw = dx2_nm2 ./ max.(vx_nm2, eps()) .+ dy2_nm2 ./ max.(vy_nm2, eps())
    A = 0.0
    B = median(z_raw) / (2 * log(2))   # start; also the scale-only refit's start if needed

    iter(A0, B0; kw...) = _iterate_calibration(reported_all, observed_all, bin_of, n_bins,
                                               vx_nm2, vy_nm2, dx2_nm2, dy2_nm2, A0, B0; kw...)
    fallback(msg, n) = _fallback_result(mean_chi2, n, n_tracks_used, frame_shifts, msg;
                                        n_tracks_filtered=n_tracks_filtered)

    fit, keep, converged = iter(A, B)
    fit isa String && return fallback(fit, count(keep))
    count(keep) < 100 &&
        return fallback("Too few pairs after chi² pair filtering ($(count(keep)) < 100)", count(keep))
    free = fit   # reported as A, B, their errors, R² and the bins

    # The free fit separates a motion term from the CRLB scale only when its slope is
    # positive and determined (B > 0, B ≥ 3 standard errors, R² ≥ 0.1). When it is not,
    # one parameter must go: either the reported precisions span too narrow a range to
    # tell the two terms apart, or a motion term dominates and the observed variance
    # barely changes with the CRLB. Both one-parameter models are fit, the CRLB scale
    # alone (σ_motion = 0) and the motion term alone (CRLB as reported, k = 1), and the
    # one that fits the bins better is applied. On a narrow range both give the right
    # total; on a wide range only the right one follows the data across the range.
    # A determined free fit with A < 0 refits the scale alone: the bound σ_motion² ≥ 0
    # binds, and keeping the slope fitted jointly with the negative intercept moves the
    # applied variance away from the data.
    determined = fit.B > 0 && fit.B >= _CAL_MIN_SLOPE_SIGNIF * fit.B_sigma &&
                 fit.r_squared >= 0.1
    note = ""
    if !determined
        sfit, skeep, sconv = iter(0.0, free.B > 0 ? free.B : B; fixA = 0.0)
        mfit, mkeep, mconv = iter(max(free.A, 0.0), 1.0; fixB = 1.0)
        sok = !(sfit isa String) && sfit.B > 0 &&
              sfit.B >= _CAL_MIN_SLOPE_SIGNIF * sfit.B_sigma
        mok = !(mfit isa String) && mfit.A > 0 &&
              mfit.A >= _CAL_MIN_SLOPE_SIGNIF * mfit.A_sigma
        if !sok && !mok
            return fallback("Neither the CRLB scale nor a motion term is determined " *
                            "(free fit: B = $(round(free.B, sigdigits=3)) ± " *
                            "$(round(free.B_sigma, sigdigits=3)), R² = " *
                            "$(round(free.r_squared, digits=3)))", count(keep))
        end
        use_scale = sok && (!mok || sfit.misfit <= mfit.misfit)
        fit, keep, conv = use_scale ? (sfit, skeep, sconv) : (mfit, mkeep, mconv)
        converged &= conv
        note = "Motion and CRLB scale not separable (free fit: B = " *
               "$(round(free.B, sigdigits=3)) ± $(round(free.B_sigma, sigdigits=3)), R² = " *
               "$(round(free.r_squared, digits=3))); applied the better one-parameter fit: " *
               (use_scale ? "CRLB scale alone, σ_motion = 0" : "motion term alone, k = 1")
        @warn "analyze_calibration: $note"
    elseif fit.A < 0
        fit, keep, conv = iter(0.0, fit.B; fixA = 0.0)
        converged &= conv
        fit isa String && return fallback(fit, count(keep))
    end
    if config.clamp_k_to_one && fit.B < 1
        fit, keep, conv = iter(max(fit.A, 0.0), 1.0; fixB = 1.0)
        converged &= conv
        fit isa String && return fallback(fit, count(keep))
        if fit.A < 0   # both bounds bind: the CRLB as reported, no motion term
            fit = merge(fit, (A = 0.0, A_sigma = 0.0))
        end
    end
    n_pairs = count(keep)
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
        n_pairs, n_tracks_used, n_tracks_filtered,
        bin_centers .* 1e-6, bin_observed .* 1e-6, copy(n_per_bin), frame_shifts,
        true, join(filter(!isempty, [note, converged ? "" :
              "Calibration fit did not converge in $(_CAL_MAX_ITER) iterations; the last fit was applied"]), "; "),
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
        new_σx_sq = ET(sigma_motion_sq + k_sq * e.σ_x^2)
        new_σy_sq = ET(sigma_motion_sq + k_sq * e.σ_y^2)
        new_σxy = ET(k_sq * e.σ_xy)

        new_emitters[i] = SMLMData.Emitter2DFit{ET}(
            e.x, e.y, e.photons, e.bg,
            sqrt(new_σx_sq), sqrt(new_σy_sq), new_σxy,
            e.σ_photons, e.σ_bg,
            e.frame, e.dataset, e.track_id, e.id
        )
    end

    return BasicSMLD(new_emitters, smld.camera, smld.n_frames,
                     smld.n_datasets, copy(smld.metadata))
end

function _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts, warning;
                          n_tracks_filtered=0)
    return CalibrationResult(
        0.0, 1.0,               # sigma_motion_nm, k_scale
        0.0, 1.0, 0.0, 0.0,    # A, B, A_sigma, B_sigma
        0.0, mean_chi2,         # r_squared, mean_chi2
        n_pairs, n_tracks_used, n_tracks_filtered,
        Float64[], Float64[], Int[],   # bin_centers, bin_observed, bin_counts
        frame_shifts,
        false, warning,
        0.0, 1.0                # A_applied, B_applied: no correction
    )
end
