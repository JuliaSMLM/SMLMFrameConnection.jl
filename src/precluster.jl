using NearestNeighbors
using Statistics
using StatsBase

"""
    precluster(smld::BasicSMLD, params::ParamStruct = ParamStruct())

Cluster localizations in `smld` based on spatiotemporal proximity thresholds.

# Description
Localizations are clustered together based on their spatiotemporal separations.
All localizations within a spatial threshold of `params.nsigmadev * mean([σ_x, σ_y])`
and a temporal threshold of `params.maxframegap` frames from one another will be
clustered together, meaning they share the same unique `track_id` value.

# Notes
Pre-clustering allows localizations observed in the same frame to be in the same
cluster. This prevents exclusion of the "correct" localization from its ideal
cluster due to a previously included "incorrect" localization.

# Returns
A new `BasicSMLD` with updated `track_id` values in the emitters.
"""
function precluster(smld::BasicSMLD{T, E}, params::ParamStruct = ParamStruct()) where {T, E<:Emitter2DFit}
    # Convert to StructArray for efficient column access
    emitters = StructArray(smld.emitters)
    n = length(emitters)

    # Sort by dataset
    sortindices_ds = sortperm(emitters.dataset)
    frames = emitters.frame[sortindices_ds]
    datasets = emitters.dataset[sortindices_ds]
    xy = hcat(emitters.x[sortindices_ds], emitters.y[sortindices_ds])'
    σ_x = emitters.σ_x[sortindices_ds]
    σ_y = emitters.σ_y[sortindices_ds]
    mean_σ = vec(mean(hcat(σ_x, σ_y), dims=2))

    # Extract parameters
    maxframegap = params.maxframegap
    nsigmadev = params.nsigmadev
    nmaxnn = params.nmaxnn

    # Initialize track IDs - each localization starts as unique cluster
    track_ids = collect(1:n)

    # Process each dataset separately
    n_per_dataset = counts(datasets)
    n_per_dataset = n_per_dataset[n_per_dataset .!= 0]
    cumulative_ds = [0; cumsum(n_per_dataset)]
    max_id = n

    for ds_idx in 1:length(unique(datasets))
        # Extract indices for current dataset
        current_idx_ds = (1:n_per_dataset[ds_idx]) .+ cumulative_ds[ds_idx]

        # Sort by frame within dataset
        sortindices_frame = sortperm(frames[current_idx_ds])
        current_idx_ds = current_idx_ds[sortindices_frame]

        frames_ds = frames[current_idx_ds]
        xy_ds = xy[:, current_idx_ds]
        mean_σ_ds = mean_σ[current_idx_ds]
        track_ids_ds = track_ids[current_idx_ds]

        # Process frames
        n_per_frame = counts(frames_ds)
        n_per_frame = n_per_frame[n_per_frame .!= 0]
        cumulative_frame = [0; cumsum(n_per_frame)]
        cluster_indices = [[idx] for idx in 1:n_per_dataset[ds_idx]]
        unique_frames = unique(frames_ds)

        for frame_idx in 1:length(unique_frames)
            current_frame = unique_frames[frame_idx]
            current_idx_frame = (1:n_per_frame[frame_idx]) .+ cumulative_frame[frame_idx]

            # Find candidate localizations for clustering (within frame gap)
            candidate_idx = findall(
                (frames_ds .>= (current_frame - maxframegap)) .&
                (frames_ds .<= current_frame)
            )

            if length(candidate_idx) < 2
                # No clustering possible - assign new unique IDs
                max_id += 1
                track_ids_ds[current_idx_frame] = (1:n_per_frame[frame_idx]) .+ max_id
                max_id += n_per_frame[frame_idx]
                continue
            end

            # Build KD-tree for nearest neighbor search
            kdtree = KDTree(xy_ds[:, candidate_idx])
            nn_indices, nn_distances = knn(
                kdtree,
                xy_ds[:, current_idx_frame],
                min(nmaxnn + 1, length(candidate_idx)),
                true
            )

            # Assign localizations to clusters based on distance cutoff
            for i in 1:n_per_frame[frame_idx]
                # Determine which candidates meet distance criterion
                σ_sum = mean_σ_ds[current_idx_frame[i]] .+ mean_σ_ds[candidate_idx[nn_indices[i]]]
                valid_nn_idx = nn_indices[i][nn_distances[i] .<= (nsigmadev * σ_sum)]

                # Update track IDs to reflect merged clusters
                update_idx = unique([
                    current_idx_frame[i];
                    candidate_idx[valid_nn_idx];
                    findall(track_ids_ds .== track_ids_ds[current_idx_frame[i]]);
                    cluster_indices[candidate_idx[valid_nn_idx]]...
                ])

                track_ids_ds[update_idx] .= minimum(track_ids_ds[update_idx])

                # Track cluster membership
                for j in update_idx
                    cluster_indices[j] = [cluster_indices[j][1]; update_idx]
                end
            end
        end

        track_ids[current_idx_ds] = track_ids_ds
    end

    # Compress track IDs to range 1:n_clusters
    compressed_track_ids = compress_connectID(track_ids)

    # Create new emitters with updated track IDs
    new_emitters = Vector{E}(undef, n)
    for i in 1:n
        idx = sortindices_ds[i]
        new_emitters[idx] = E(
            emitters.x[idx],
            emitters.y[idx],
            emitters.photons[idx],
            emitters.bg[idx],
            emitters.σ_x[idx],
            emitters.σ_y[idx],
            emitters.σ_photons[idx],
            emitters.σ_bg[idx],
            emitters.frame[idx],
            emitters.dataset[idx],
            compressed_track_ids[i],
            emitters.id[idx]
        )
    end

    return BasicSMLD(new_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end
