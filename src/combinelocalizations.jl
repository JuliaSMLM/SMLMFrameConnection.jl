using SMLMData
using StatsBase

"""
    smld_combined = combinelocalizations(smld::BasicSMLD{T,E};
        track_length::Union{Tuple{Float64,Float64},Nothing}=nothing) where {T, E<:SMLMData.AbstractEmitter}

Combine clustered localizations in `smld` into higher precision localizations.

# Description
This function combines localizations in `smld` that share the same value of
track_id.  Localizations are combined assuming they arose from
independent measurements of the same position with Gaussian errors.

`track_length` is an inclusive `(min, max)` range on the number of localizations
per track: a track is emitted only if `min <= nperID <= max`. `nothing` (the
default) disables filtering. This counts localizations sharing a `track_id`, not
temporal frame-span. Use `(2.0, Inf)` to drop single-frame blinks, or a finite
upper bound like `(2.0, 50.0)` to also drop over-long tracks (e.g. fiducials or
sticky docking strands in dense DNA-PAINT data).
"""
function combinelocalizations(smld::BasicSMLD{T,E};
    track_length::Union{Tuple{Float64,Float64},Nothing}=nothing) where {T, E<:SMLMData.AbstractEmitter}
    # Extract arrays from emitters
    emitters = smld.emitters
    # Nothing to combine: return the (empty) input unchanged. The reductions
    # below (`counts`, `first(emitters)`) assume at least one localization.
    isempty(emitters) && return smld
    connectID = [e.track_id for e in emitters]
    sortindices = sortperm(connectID)
    connectID = connectID[sortindices]
    x = [e.x for e in emitters][sortindices]
    y = [e.y for e in emitters][sortindices]
    σ_x = [e.σ_x for e in emitters][sortindices]
    σ_y = [e.σ_y for e in emitters][sortindices]
    σ_xy = [e.σ_xy for e in emitters][sortindices]
    photons = [e.photons for e in emitters][sortindices]
    σ_photons = [e.σ_photons for e in emitters][sortindices]
    bg = [e.bg for e in emitters][sortindices]
    σ_bg = [e.σ_bg for e in emitters][sortindices]
    framenum = [e.frame for e in emitters][sortindices]
    datasetnum = [e.dataset for e in emitters][sortindices]
    id = [e.id for e in emitters][sortindices]

    # Loop over clusters and combine their localization data as appropriate.
    nperID = counts(connectID)
    nperID = nperID[nperID .!= 0]
    # ncumulative is built over the FULL set of present clusters so its offsets
    # span every localization in the sorted arrays. The length filter only
    # selects which clusters to emit (`keep`); it must not alter these offsets.
    ncumulative = [0; cumsum(nperID)]
    keep = if track_length === nothing
        collect(1:length(nperID))
    else
        lo, hi = track_length
        findall(n -> lo <= n <= hi, nperID)
    end
    nclusters = length(keep)

    # Pre-allocate arrays for combined values
    # Use emitter's native precision, not SMLD's type parameter
    ET = typeof(first(emitters).x)
    x_combined = Vector{ET}(undef, nclusters)
    y_combined = Vector{ET}(undef, nclusters)
    σ_x_combined = Vector{ET}(undef, nclusters)
    σ_y_combined = Vector{ET}(undef, nclusters)
    σ_xy_combined = Vector{ET}(undef, nclusters)
    photons_combined = Vector{ET}(undef, nclusters)
    σ_photons_combined = Vector{ET}(undef, nclusters)
    bg_combined = Vector{ET}(undef, nclusters)
    σ_bg_combined = Vector{ET}(undef, nclusters)
    connectID_combined = Vector{Int}(undef, nclusters)
    framenum_combined = Vector{Int}(undef, nclusters)
    datasetnum_combined = Vector{Int}(undef, nclusters)
    id_combined = Vector{Int}(undef, nclusters)

    for (out_idx, nn) in enumerate(keep)
        indices = (1:nperID[nn]) .+ ncumulative[nn]

        # Combine positions using precision-weighted mean with full 2x2 covariance.
        # For each measurement i: Σᵢ = [σ_x² σ_xy; σ_xy σ_y²], Pᵢ = Σᵢ⁻¹
        # Combined: Σ_comb = (Σ Pᵢ)⁻¹, pos_comb = Σ_comb * Σ(Pᵢ * posᵢ)
        # Accumulate in Float64 to avoid precision loss (eps(Float32)=1.2e-7
        # is larger than typical determinants for SMLM sigmas ~0.006 μm)
        P_sum = zeros(Float64, 2, 2)
        μ_weighted = zeros(Float64, 2)
        for i in indices
            sx2 = Float64(σ_x[i])^2
            sy2 = Float64(σ_y[i])^2
            sxy = Float64(σ_xy[i])
            det_i = sx2 * sy2 - sxy^2
            det_i = max(det_i, eps(Float64))  # Ensure positive definite
            # Precision matrix: P = [σ_y² -σ_xy; -σ_xy σ_x²] / det
            P11 = sy2 / det_i
            P22 = sx2 / det_i
            P12 = -sxy / det_i
            P_sum[1,1] += P11
            P_sum[2,2] += P22
            P_sum[1,2] += P12
            P_sum[2,1] += P12
            μ_weighted[1] += P11 * Float64(x[i]) + P12 * Float64(y[i])
            μ_weighted[2] += P12 * Float64(x[i]) + P22 * Float64(y[i])
        end

        # Combined covariance = inverse of summed precision
        det_P = P_sum[1,1] * P_sum[2,2] - P_sum[1,2]^2
        Σ_comb_11 = P_sum[2,2] / det_P  # σ_x²
        Σ_comb_22 = P_sum[1,1] / det_P  # σ_y²
        Σ_comb_12 = -P_sum[1,2] / det_P  # σ_xy

        # Combined position (cast back to emitter precision)
        x_combined[out_idx] = ET(Σ_comb_11 * μ_weighted[1] + Σ_comb_12 * μ_weighted[2])
        y_combined[out_idx] = ET(Σ_comb_12 * μ_weighted[1] + Σ_comb_22 * μ_weighted[2])
        σ_x_combined[out_idx] = ET(sqrt(Σ_comb_11))
        σ_y_combined[out_idx] = ET(sqrt(Σ_comb_22))
        σ_xy_combined[out_idx] = ET(Σ_comb_12)

        # Photons and background: sum (independent measurements)
        photons_combined[out_idx] = sum(photons[indices])
        σ_photons_combined[out_idx] = sqrt(sum(σ_photons[indices] .^ 2))
        bg_combined[out_idx] = sum(bg[indices])
        σ_bg_combined[out_idx] = sqrt(sum(σ_bg[indices] .^ 2))

        connectID_combined[out_idx] = connectID[indices[1]]
        framenum_combined[out_idx] = framenum[indices[1]]
        datasetnum_combined[out_idx] = datasetnum[indices[1]]
        id_combined[out_idx] = id[indices[1]]
    end

    # Create new emitters for combined results
    combined_emitters = Vector{SMLMData.Emitter2DFit{ET}}(undef, nclusters)
    for nn in 1:nclusters
        combined_emitters[nn] = SMLMData.Emitter2DFit{ET}(
            x_combined[nn], y_combined[nn],
            photons_combined[nn], bg_combined[nn],
            σ_x_combined[nn], σ_y_combined[nn],
            σ_xy_combined[nn],
            σ_photons_combined[nn], σ_bg_combined[nn],
            framenum_combined[nn], datasetnum_combined[nn],
            connectID_combined[nn], id_combined[nn]
        )
    end

    smld_combined = BasicSMLD(combined_emitters, smld.camera, smld.n_frames,
                              smld.n_datasets, copy(smld.metadata))

    return smld_combined
end
