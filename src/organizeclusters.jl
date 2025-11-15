"""
    organizeclusters(smld::BasicSMLD)

Organize localizations by cluster into a vector of matrices.

# Description
Localizations in `smld`, organized via their `track_id` values, are organized
into a vector of matrices `clusterdata`. Each matrix contains the data for one
cluster, with columns: [x, y, σ_x, σ_y, frame, dataset, track_id, original_index].

The 8th column (original_index) stores the indices into the original emitter array,
used for mapping LAP solutions back to emitters.

# Returns
A `Vector{Matrix{Float32}}` where each matrix represents one cluster's localization data.
"""
function organizeclusters(smld::BasicSMLD{T, E}) where {T, E<:Emitter2DFit}
    # Convert to StructArray for efficient column access
    emitters = StructArray(smld.emitters)

    # Sort by track_id
    track_ids = emitters.track_id
    sortindices = sortperm(track_ids)

    # Organize into combined data matrix
    # Column 8 stores original indices for mapping back after LAP
    combined_data = hcat(
        emitters.x[sortindices],
        emitters.y[sortindices],
        emitters.σ_x[sortindices],
        emitters.σ_y[sortindices],
        Float64.(emitters.frame[sortindices]),
        Float64.(emitters.dataset[sortindices]),
        Float64.(track_ids[sortindices]),
        Float64.(sortindices)
    )

    # Get cluster counts
    track_ids_sorted = track_ids[sortindices]
    n_per_cluster = counts(track_ids_sorted)
    n_per_cluster = n_per_cluster[n_per_cluster .!= 0]
    cumulative = [0; cumsum(n_per_cluster)]
    n_clusters = length(n_per_cluster)

    # Extract data for each cluster
    clusterdata = Vector{Matrix{Float32}}(undef, n_clusters)
    for cluster_idx in 1:n_clusters
        cluster_indices = (1:n_per_cluster[cluster_idx]) .+ cumulative[cluster_idx]
        clusterdata[cluster_idx] = Float32.(combined_data[cluster_indices, :])
    end

    return clusterdata
end
