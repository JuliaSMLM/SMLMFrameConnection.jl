"""
    connect1DS!(smld::BasicSMLD, dataset::Int; maxframegap::Int = 5)

Define the "ideal" frame-connection result for simulated `smld` with one dataset.

# Description
This is a helper function which defines the ideal FC result for a single dataset.
The user should probably not bother using this, and instead should call `defineidealFC()`.

Localizations belonging to the same emitter (as indicated by their `.id` field) are
assigned the same `track_id` unless they are separated by more than `maxframegap` frames.

# Arguments
- `smld`: BasicSMLD with emitter `.id` values populated to indicate true emitter membership
- `dataset`: Dataset number to be connected
- `maxframegap`: Maximum frame gap allowed between localizations connected in the "ideal" result

# Returns
Modified `smld` with updated `track_id` values for the specified dataset.
"""
function connect1DS!(smld::BasicSMLD{T, E}, dataset::Int; maxframegap::Int = 5) where {T, E<:Emitter2DFit}
    # Get maximum existing track_id
    max_track_id = maximum(e.track_id for e in smld.emitters; init=0)

    # Create new emitters array
    new_emitters = Vector{E}(undef, length(smld.emitters))

    # Process emitters
    for (idx, emitter) in enumerate(smld.emitters)
        if emitter.dataset != dataset
            # Keep track_id unchanged for other datasets
            new_emitters[idx] = emitter
        else
            # Will be updated below
            new_emitters[idx] = set_track_id(emitter, emitter.track_id + max_track_id)
        end
    end

    # Group by emitter ID within this dataset
    dataset_indices = findall(e -> e.dataset == dataset, new_emitters)
    emitter_ids = unique(new_emitters[i].id for i in dataset_indices)

    for emitter_id in emitter_ids
        # Find all localizations of this emitter in current dataset
        current_indices = findall(i -> new_emitters[i].dataset == dataset &&
                                       new_emitters[i].id == emitter_id,
                                  dataset_indices)
        current_indices = dataset_indices[current_indices]

        if isempty(current_indices)
            continue
        end

        # Sort by frame
        frames = [new_emitters[i].frame for i in current_indices]
        sort_idx = sortperm(frames)
        sorted_indices = current_indices[sort_idx]
        sorted_frames = frames[sort_idx]

        # Check if all within maxframegap
        frame_diffs = diff(sorted_frames)
        if all(frame_diffs .<= maxframegap)
            # All localizations can share the same track_id (already set)
            continue
        end

        # Split into separate blinking events based on frame gaps
        for i in 1:length(frame_diffs)
            if frame_diffs[i] <= maxframegap
                # Connect these localizations
                new_emitters[sorted_indices[i+1]] = set_track_id(
                    new_emitters[sorted_indices[i+1]],
                    new_emitters[sorted_indices[i]].track_id
                )
            else
                # Start new blinking event
                max_track_id += 1
                new_emitters[sorted_indices[i+1]] = set_track_id(
                    new_emitters[sorted_indices[i+1]],
                    max_track_id
                )
            end
        end
    end

    return BasicSMLD(new_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""
    connect1DS(smld::BasicSMLD, dataset::Int; maxframegap::Int = 5)

Non-mutating version of `connect1DS!`.
"""
function connect1DS(smld::BasicSMLD{T, E}, dataset::Int; maxframegap::Int = 5) where {T, E<:Emitter2DFit}
    return connect1DS!(deepcopy(smld), dataset; maxframegap=maxframegap)
end

"""
    defineidealFC!(smld::BasicSMLD; maxframegap::Int = 5)

Define the "ideal" frame-connection result for a simulated `smld`.

# Description
For simulated BasicSMLD structures where the emitter `.id` field indicates true
emitter membership, this function generates an "ideal" FC result which combines
all blinking events that appeared with frame gaps less than `maxframegap`.

Note: For very high duty cycles, multiple blinking events might be mistakenly
combined if the emitter blinks back on within `maxframegap` frames. Localizations
are not allowed to be connected across datasets.

# Arguments
- `smld`: BasicSMLD with emitter `.id` values populated to indicate true emitter membership
- `maxframegap`: Maximum frame gap allowed between connected localizations (default = 5)

# Returns
- `smld`: Input SMLD with `track_id` modified to reflect ideal blinking event connections
- `smld_combined`: Result of combining localizations with same `track_id`
"""
function defineidealFC!(smld::BasicSMLD{T, E}; maxframegap::Int = 5) where {T, E<:Emitter2DFit}
    # Process each dataset
    for dataset_id in unique(e.dataset for e in smld.emitters)
        smld = connect1DS!(smld, dataset_id; maxframegap=maxframegap)
    end

    # Compress track IDs to sequential 1:n_clusters
    track_ids = [e.track_id for e in smld.emitters]
    compress_connectID!(track_ids)

    # Update emitters with compressed track IDs
    new_emitters = [set_track_id(smld.emitters[i], track_ids[i]) for i in 1:length(smld.emitters)]
    smld = BasicSMLD(new_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)

    # Combine localizations
    smld_combined = combinelocalizations(smld)

    return smld, smld_combined
end

"""
    defineidealFC(smld::BasicSMLD; maxframegap::Int = 5)

Non-mutating version of `defineidealFC!`.
"""
function defineidealFC(smld::BasicSMLD{T, E}; maxframegap::Int = 5) where {T, E<:Emitter2DFit}
    return defineidealFC!(deepcopy(smld); maxframegap=maxframegap)
end

"""
    set_track_id(emitter::Emitter2DFit, track_id::Int)

Create a new emitter with updated track_id.
"""
function set_track_id(emitter::Emitter2DFit{T}, track_id::Int) where T
    return Emitter2DFit{T}(
        emitter.x,
        emitter.y,
        emitter.photons,
        emitter.bg,
        emitter.σ_x,
        emitter.σ_y,
        emitter.σ_photons,
        emitter.σ_bg,
        emitter.frame,
        emitter.dataset,
        track_id,
        emitter.id
    )
end
