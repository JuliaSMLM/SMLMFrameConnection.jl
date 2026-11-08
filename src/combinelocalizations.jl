using SMLMData
using StatsBase

"""
    smld_combined = combinelocalizations(smld::BasicSMLD{T,E};
        track_length::Union{Tuple{Float64,Float64},Nothing}=nothing) where {T, E<:SMLMData.AbstractEmitter}

Combine clustered localizations in `smld` into higher precision localizations.

# Description
This function combines localizations in `smld` that share the same value of
track_id within one dataset.  Localizations are combined assuming they arose from
independent measurements of the same position with Gaussian errors. A localization
with `track_id == 0` (SMLMData's "unlinked") is passed through on its own. Each
combined localization takes frame, dataset and id from its track's earliest
localization; the result does not depend on the order of the input.

`track_length` is an inclusive `(min, max)` range on the number of localizations
per track: a track is emitted only if `min <= nperID <= max`. `nothing` (the
default) disables filtering. This counts localizations sharing a `track_id`, not
temporal frame-span. Use `(2.0, Inf)` to drop single-frame blinks, or a finite
upper bound like `(2.0, 50.0)` to also drop over-long tracks (e.g. fiducials or
sticky docking strands in dense DNA-PAINT data).
"""
function combinelocalizations(smld::BasicSMLD{T,E};
    track_length::Union{Tuple{Float64,Float64},Nothing}=nothing) where {T, E<:SMLMData.AbstractEmitter}
    _check_2d(smld, "combinelocalizations")
    _check_inputs(smld, "combinelocalizations"; frames = false)
    # Extract arrays from emitters
    emitters = smld.emitters
    # Nothing to combine: return the (empty) input unchanged. The reductions
    # below (`first(emitters)`) assume at least one localization.
    isempty(emitters) && return smld
    connectID = [e.track_id for e in emitters]
    # Order by (track, dataset, frame, x, y), ties broken by the other fields, so the
    # result does not depend on the input order: each combined localization takes frame,
    # dataset and id from its track's first localization in time.
    sortindices = sortperm(emitters; by = e -> (e.track_id, _order_key(e)))
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

    # Groups: runs of one track_id within one dataset; every track_id-0 localization
    # is a group of its own (unlinked).
    nperID = Int[]
    for i in eachindex(connectID)
        newgroup = i == 1 || connectID[i] == 0 || connectID[i] != connectID[i-1] ||
                   datasetnum[i] != datasetnum[i-1]
        if newgroup
            push!(nperID, 1)
        else
            nperID[end] += 1
        end
    end
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

        # A group of one (an unlinked localization, or a one-localization track) is its
        # own combination: copied, not round-tripped through the precision matrix.
        if nperID[nn] == 1
            i = indices[1]
            x_combined[out_idx], y_combined[out_idx] = x[i], y[i]
            σ_x_combined[out_idx], σ_y_combined[out_idx] = σ_x[i], σ_y[i]
            σ_xy_combined[out_idx] = σ_xy[i]
            photons_combined[out_idx], σ_photons_combined[out_idx] = photons[i], σ_photons[i]
            bg_combined[out_idx], σ_bg_combined[out_idx] = bg[i], σ_bg[i]
            connectID_combined[out_idx] = connectID[i]
            framenum_combined[out_idx] = framenum[i]
            datasetnum_combined[out_idx] = datasetnum[i]
            id_combined[out_idx] = id[i]
            continue
        end

        # Combine positions using precision-weighted mean with full 2x2 covariance.
        # For each measurement i: Σᵢ = [σ_x² σ_xy; σ_xy σ_y²], Pᵢ = Σᵢ⁻¹
        # Combined: Σ_comb = (Σ Pᵢ)⁻¹, pos_comb = Σ_comb * Σ(Pᵢ * posᵢ)
        # In Float64, with the summed precision's determinant carried accurately
        # (_precision_mean); cast back to emitter precision with the covariance kept
        # positive definite through the rounding.
        μx, μy, Σ11, Σ22, Σ12 = _precision_mean(x, y, σ_x, σ_y, σ_xy, indices)
        x_combined[out_idx] = ET(μx)
        y_combined[out_idx] = ET(μy)
        σ_x_combined[out_idx], σ_y_combined[out_idx], σ_xy_combined[out_idx] =
            _round_cov(ET, sqrt(Σ11), sqrt(Σ22), Σ12)

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
