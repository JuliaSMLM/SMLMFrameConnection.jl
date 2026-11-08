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
# Returns a NamedTuple, or a String naming why no fit was possible.
function _binned_wls(reported, observed, bin_of, keep, n_bins, A, B)
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

    X = hcat(ones(n_valid_bins), bin_centers)
    W = Diagonal(weights)
    XtWX = X' * W * X
    XtWy = X' * W * bin_observed
    det_XtWX = XtWX[1,1] * XtWX[2,2] - XtWX[1,2]^2
    abs(det_XtWX) < eps() && return "Singular WLS matrix"

    coeffs = XtWX \ XtWy
    y_pred = X * coeffs
    residuals = bin_observed .- y_pred
    ss_res = sum(weights .* residuals .^ 2)
    ss_tot = sum(weights .* (bin_observed .- mean(bin_observed)) .^ 2)
    r_squared = ss_tot > 0 ? 1.0 - ss_res / ss_tot : NaN
    mse = ss_res / max(1, n_valid_bins - 2)
    var_coeffs = mse * inv(XtWX)

    return (A = coeffs[1], B = coeffs[2],
            A_sigma = sqrt(max(0.0, var_coeffs[1, 1])),
            B_sigma = sqrt(max(0.0, var_coeffs[2, 2])),
            r_squared = r_squared, bin_centers = bin_centers,
            bin_observed = bin_observed, n_per_bin = n_per_bin)
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

Falls back gracefully (returns `calibration_applied=false`) if fewer than 100
frame-to-frame pairs are available, if the fit quality is poor (R² < 0.1), or if the
slope `B` is not positive and at least 3 standard errors (the reported precisions span
too narrow a range to separate the CRLB scale from the motion term).
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
    B = median(z_raw) / (2 * log(2))

    fit = "No calibration fit was attempted"
    keep = trues(length(dx_vec))
    converged = false
    for _ in 1:_CAL_MAX_ITER
        sx = max.(A .+ B .* vx_nm2, eps())
        sy = max.(A .+ B .* vy_nm2, eps())
        keep = (dx2_nm2 ./ sx .+ dy2_nm2 ./ sy) .<= _CAL_CHI2_CUT
        fit = _binned_wls(reported_all, observed_all, bin_of, keep, n_bins, A, B)
        fit isa String && break
        done = abs(fit.B - B) <= 1e-6 * max(1.0, abs(B)) &&
               abs(fit.A - A) <= 1e-6 * max(1.0, mean(reported_all))
        A, B = fit.A, fit.B
        converged = done
        done && break
    end
    n_pairs = count(keep)

    if n_pairs < 100
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Too few pairs after chi² pair filtering ($n_pairs < 100)";
                                n_tracks_filtered=n_tracks_filtered)
    end
    if fit isa String
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts, fit;
                                n_tracks_filtered=n_tracks_filtered)
    end

    (; bin_centers, bin_observed, n_per_bin, r_squared, A_sigma, B_sigma) = fit
    if r_squared < 0.1
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Poor fit quality (R² = $(round(r_squared, digits=3)))";
                                n_tracks_filtered=n_tracks_filtered)
    end
    # The slope must be positive and determined: when the reported precisions span a
    # narrow range, intercept and slope trade off and the fit can land on any slope,
    # including a negative one with a large intercept (k = 0 and a large motion term).
    if !(B > 0 && B >= _CAL_MIN_SLOPE_SIGNIF * B_sigma)
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "CRLB scale not determined: slope B = $(round(B, sigdigits=3)) " *
                                "± $(round(B_sigma, sigdigits=3)) (the reported precisions " *
                                "may span too narrow a range)";
                                n_tracks_filtered=n_tracks_filtered)
    end

    # Extract parameters -- fit was in nm², convert to μm² for storage.
    # A is the pair intercept 2·σ_motion² (nm²), B is k² (dimensionless).
    sigma_motion_nm = sqrt(max(0.0, A) / 2)
    k_sq = max(0.0, B)
    k_scale = sqrt(k_sq)
    if config.clamp_k_to_one
        k_scale = max(k_scale, 1.0)
    end

    # Convert A, A_sigma, bin data from nm² to μm² for CalibrationResult
    A_um2 = A * 1e-6
    A_sigma_um2 = A_sigma * 1e-6
    bin_centers_um2 = bin_centers .* 1e-6
    bin_observed_um2 = bin_observed .* 1e-6

    return CalibrationResult(
        sigma_motion_nm, k_scale,
        A_um2, B, A_sigma_um2, B_sigma, r_squared, mean_chi2,
        n_pairs, n_tracks_used, n_tracks_filtered,
        bin_centers_um2, bin_observed_um2, copy(n_per_bin), frame_shifts,
        true, converged ? "" :
              "Calibration fit did not converge in $(_CAL_MAX_ITER) iterations; the last fit was applied"
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
        false, warning
    )
end
