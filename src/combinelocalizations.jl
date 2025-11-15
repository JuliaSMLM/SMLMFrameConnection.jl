using StatsBase

"""
    combinelocalizations(smld::BasicSMLD)

Combine clustered localizations in `smld` into higher precision localizations.

# Description
This function combines localizations in `smld` that share the same `track_id`.
Localizations are combined assuming they arose from independent measurements of
the same position with Gaussian errors, using inverse variance weighting.

# Returns
A new `BasicSMLD` with one emitter per unique `track_id`, containing combined
localization data.
"""
function combinelocalizations(smld::BasicSMLD{T, E}) where {T, E<:Emitter2DFit}
    # Convert to StructArray for efficient column access
    emitters = StructArray(smld.emitters)

    # Sort by track_id
    sortindices = sortperm(emitters.track_id)

    # Get counts for each track_id
    track_ids_sorted = emitters.track_id[sortindices]
    n_per_track = counts(track_ids_sorted)
    n_per_track = n_per_track[n_per_track .!= 0]
    cumulative = [0; cumsum(n_per_track)]
    n_clusters = length(n_per_track)

    # Combine emitters for each cluster
    combined_emitters = Vector{E}(undef, n_clusters)

    for cluster_idx in 1:n_clusters
        indices = sortindices[(1:n_per_track[cluster_idx]) .+ cumulative[cluster_idx]]
        combined_emitters[cluster_idx] = combine_emitter_group(emitters, indices, E)
    end

    return BasicSMLD(combined_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""
    combine_emitter_group(emitters::StructArray, indices::Vector{Int}, ::Type{E})

Combine a group of emitters using inverse variance weighted averaging.

Positions are weighted by inverse variance (1/σ²). Photons and backgrounds are
summed. Uncertainties are propagated according to error propagation rules.
"""
function combine_emitter_group(emitters::StructArray, indices::Vector{Int}, ::Type{E}) where {E<:Emitter2DFit}
    # Inverse variance weights for position
    w_x = 1.0 ./ emitters.σ_x[indices].^2
    w_y = 1.0 ./ emitters.σ_y[indices].^2

    # Weighted mean positions
    x_combined = sum(emitters.x[indices] .* w_x) / sum(w_x)
    y_combined = sum(emitters.y[indices] .* w_y) / sum(w_y)

    # Combined uncertainties (inverse variance weighting)
    σ_x_combined = sqrt(1.0 / sum(w_x))
    σ_y_combined = sqrt(1.0 / sum(w_y))

    # Sum photons and backgrounds
    photons_combined = sum(emitters.photons[indices])
    bg_combined = sum(emitters.bg[indices])

    # Propagate uncertainties (uncorrelated errors)
    σ_photons_combined = sqrt(sum(emitters.σ_photons[indices].^2))
    σ_bg_combined = sqrt(sum(emitters.σ_bg[indices].^2))

    # Use metadata from first emitter in cluster
    return E(
        x_combined,
        y_combined,
        photons_combined,
        bg_combined,
        σ_x_combined,
        σ_y_combined,
        σ_photons_combined,
        σ_bg_combined,
        emitters.frame[indices[1]],
        emitters.dataset[indices[1]],
        emitters.track_id[indices[1]],
        0  # Reset ID for combined emitter
    )
end
