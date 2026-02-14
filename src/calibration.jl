# Uncertainty calibration for frame-connected tracks.
#
# Analyzes frame-to-frame jitter within connected tracks to estimate motion
# variance (σ_motion²) and CRLB scale factor (k²), then applies corrected
# uncertainties: Σ_corrected = σ_motion² I + k² Σ_CRLB

"""
    analyze_calibration(smld::BasicSMLD, config::CalibrationConfig) -> CalibrationResult

Analyze frame-to-frame jitter in connected tracks to estimate uncertainty
calibration parameters.

Collects consecutive-frame pairs within each track, computes observed positional
variance vs CRLB variance, and fits the model: `observed_var = A + B * CRLB_var`
using weighted least squares on binned data.

Falls back gracefully (returns `calibration_applied=false`) if fewer than 100
frame-to-frame pairs are available or if the fit quality is poor (R² < 0.1).
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

    # Per-pair chi² filtering: remove pairs with chi² > 6.0 (double-fits, mismatches)
    # This is critical -- reduces scatter in bins before fitting
    chi2_per_pair = dx_vec .^ 2 ./ max.(var_x_vec, eps()) .+
                    dy_vec .^ 2 ./ max.(var_y_vec, eps())
    chi2_threshold = 6.0  # ~99.7% of chi²(2) is below this
    good_mask = chi2_per_pair .<= chi2_threshold

    n_pairs_filtered = sum(.!good_mask)
    dx_vec = dx_vec[good_mask]
    dy_vec = dy_vec[good_mask]
    var_x_vec = var_x_vec[good_mask]
    var_y_vec = var_y_vec[good_mask]
    n_pairs = length(dx_vec)

    if n_pairs < 100
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Too few pairs after chi² pair filtering ($n_pairs < 100)";
                                n_tracks_filtered=n_tracks_filtered)
    end

    # Per-axis averaging and convert to nm² for fitting
    reported_var_nm2 = ((var_x_vec .+ var_y_vec) ./ 2) .* 1e6
    observed_var_nm2 = ((dx_vec .^ 2 .+ dy_vec .^ 2) ./ 2) .* 1e6

    # Equal-width bins with minimum 5 points per bin
    x_min, x_max = extrema(reported_var_nm2)
    x_range = x_max - x_min
    n_bins = 20
    bin_width = x_range / n_bins

    bin_centers = Float64[]
    bin_observed = Float64[]
    bin_var = Float64[]
    n_per_bin = Int[]

    for i in 1:n_bins
        bin_lo = x_min + (i - 1) * bin_width
        bin_hi = x_min + i * bin_width

        if i == n_bins
            mask = (reported_var_nm2 .>= bin_lo) .& (reported_var_nm2 .<= bin_hi)
        else
            mask = (reported_var_nm2 .>= bin_lo) .& (reported_var_nm2 .< bin_hi)
        end

        bin_x = reported_var_nm2[mask]
        bin_y = observed_var_nm2[mask]

        if length(bin_x) >= 5
            push!(bin_centers, mean(bin_x))
            push!(bin_observed, mean(bin_y))
            push!(bin_var, var(bin_y))
            push!(n_per_bin, length(bin_x))
        end
    end

    n_valid_bins = length(bin_centers)
    if n_valid_bins < 3
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Too few valid bins ($n_valid_bins < 3)";
                                n_tracks_filtered=n_tracks_filtered)
    end

    # WLS on binned data: weight = 1/var_of_bin_mean, normalized
    var_of_bin_mean = bin_var ./ n_per_bin
    weights = [v > 0 ? 1.0 / v : 0.0 for v in var_of_bin_mean]

    if sum(weights) ≈ 0
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "All bin weights are zero";
                                n_tracks_filtered=n_tracks_filtered)
    end

    weights = weights ./ sum(weights) .* n_valid_bins

    X = hcat(ones(n_valid_bins), bin_centers)
    W = Diagonal(weights)
    XtWX = X' * W * X
    XtWy = X' * W * bin_observed

    det_XtWX = XtWX[1,1] * XtWX[2,2] - XtWX[1,2]^2
    if abs(det_XtWX) < eps()
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Singular WLS matrix";
                                n_tracks_filtered=n_tracks_filtered)
    end

    coeffs = XtWX \ XtWy
    A = coeffs[1]  # σ_motion² in nm²
    B = coeffs[2]  # k²

    # R² and standard errors
    y_pred = X * coeffs
    residuals = bin_observed .- y_pred
    ss_res = sum(weights .* residuals .^ 2)
    ss_tot = sum(weights .* (bin_observed .- mean(bin_observed)) .^ 2)
    r_squared = ss_tot > 0 ? 1.0 - ss_res / ss_tot : NaN

    if r_squared < 0.1
        return _fallback_result(mean_chi2, n_pairs, n_tracks_used, frame_shifts,
                                "Poor fit quality (R² = $(round(r_squared, digits=3)))";
                                n_tracks_filtered=n_tracks_filtered)
    end

    mse = ss_res / max(1, n_valid_bins - 2)
    var_coeffs = mse * inv(XtWX)
    A_sigma = sqrt(max(0.0, var_coeffs[1, 1]))
    B_sigma = sqrt(max(0.0, var_coeffs[2, 2]))

    # Extract parameters -- fit was in nm², convert to μm² for storage
    # A is σ_motion² (nm²), B is k² (dimensionless)
    sigma_motion_nm = sqrt(max(0.0, A))
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
        bin_centers_um2, bin_observed_um2, frame_shifts,
        true, ""
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
        Float64[], Float64[],   # bin_centers, bin_observed
        frame_shifts,
        false, warning
    )
end
