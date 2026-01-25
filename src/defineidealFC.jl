using SMLMData

"""
    connect1DS(smld::BasicSMLD{T,E}, dataset::Int,
               connectID::Vector{Int}, maxID::Int;
               maxframegap::Int = 5) where {T, E<:SMLMData.AbstractEmitter}

Define the "ideal" frame-connection result for simulated `smld` with one dataset.

# Description
This is a helper function which defines the ideal FC result for a single
dataset.  The user should probably not bother using this, and instead should
call defineidealFC().

# Inputs
- `smld`: BasicSMLD with track_id populated to indicate emitter ID.
- `dataset`: Dataset number to be connected.
- `connectID`: Current connection IDs for all emitters.
- `maxID`: Current maximum ID value.
- `maxframegap`: Maximum frame gap allowed between localizations connected in
                 the "ideal" result.

# Returns
- Updated connectID array and new maxID.
"""
function connect1DS(smld::BasicSMLD{T,E}, dataset::Int,
                    connectID::Vector{Int}, maxID::Int;
                    maxframegap::Int = 5) where {T, E<:SMLMData.AbstractEmitter}
    emitters = smld.emitters

    # Find indices for current dataset
    currentDS = findall(e -> e.dataset == dataset, emitters)

    if isempty(currentDS)
        return connectID, maxID
    end

    # Extract relevant data for current dataset
    connectID_DS = connectID[currentDS] .+ maxID  # offset to ensure uniqueness
    maxID = maximum(connectID_DS)  # Update maxID to avoid collision with new IDs
    framenum_DS = [emitters[i].frame for i in currentDS]

    # Loop through associated localizations and combine them as appropriate.
    for ii in unique(connectID_DS)
        # Determine which localizations belong to the ii-th emitter.
        local_inds = findall(connectID_DS .== ii)

        if isempty(local_inds)
            continue
        end

        # Sort these localizations with respect to their frame.
        frames = framenum_DS[local_inds]
        sortinds = sortperm(frames)
        sortedframes = frames[sortinds]
        sorted_local_inds = local_inds[sortinds]

        # If all of these localizations are within the framegap, no action is
        # needed (they already share the same connectID).
        framediff = diff(sortedframes)
        if all(framediff .<= maxframegap)
            continue
        end

        # Determine which localizations we can combine.
        for ff = 1:length(framediff)
            if framediff[ff] <= maxframegap
                # Connect these localizations.
                connectID_DS[sorted_local_inds[ff+1]] =
                    connectID_DS[sorted_local_inds[ff]]
            else
                # Localization ff+1 should be a new blinking event.
                maxID += 1
                connectID_DS[sorted_local_inds[ff+1]] = maxID
            end
        end
    end

    # Update connectID for the current dataset.
    connectID[currentDS] = connectID_DS

    return connectID, maxID
end

"""
    smld_connected, smld_combined = defineidealFC(
        smld::BasicSMLD{T,E};
        maxframegap::Int = 5) where {T, E<:SMLMData.AbstractEmitter}

Define the "ideal" frame-connection result for a simulated `smld`.

# Description
This function defines the "ideal" frame connection result from a simulation.
That is to say, for a simulated BasicSMLD structure `smld` with `track_id` field
populated to indicate emitter membership of localizations, this function will
generate an "ideal" FC result which combines all blinking events that appeared
with frame gaps less than `maxframegap` of one another.  Note that for very
high duty cycles, multiple blinking events might be mistakingly combined by
this method (i.e., if the emitter blinks back on within `maxframegap` frames
of its previous blink).  Note that localizations are not allowed to be
connected across datasets.

# Inputs
- `smld`: BasicSMLD with track_id populated to indicate emitter ID.
- `maxframegap`: Maximum frame gap allowed between localizations connected in
                 the "ideal" result.

# Outputs
- `smld_connected`: Copy of the input `smld` with track_id modified to reflect
                    blinking event ID.
- `smld_combined`: Ideal frame-connection result with localizations combined.
"""
function defineidealFC(smld::BasicSMLD{T,E};
                       maxframegap::Int = 5) where {T, E<:SMLMData.AbstractEmitter}
    emitters = smld.emitters

    # Initialize connectID from track_id
    connectID = [e.track_id for e in emitters]
    maxID = maximum(connectID)

    # Get unique datasets
    datasets = unique(e.dataset for e in emitters)

    # Loop through datasets and combine localizations as appropriate.
    for ds in datasets
        connectID, maxID = connect1DS(smld, ds, connectID, maxID; maxframegap = maxframegap)
    end

    # Compress connectID
    compress_connectID!(connectID)

    # Create new emitters with updated track_id
    # Use emitter's native precision, not SMLD's type parameter
    ET = typeof(first(emitters).x)
    new_emitters = Vector{SMLMData.Emitter2DFit{ET}}(undef, length(emitters))
    for i in 1:length(emitters)
        e = emitters[i]
        new_emitters[i] = SMLMData.Emitter2DFit{ET}(
            e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg,
            e.frame, e.dataset, connectID[i], e.id
        )
    end

    smld_connected = BasicSMLD(new_emitters, smld.camera, smld.n_frames,
                                smld.n_datasets, copy(smld.metadata))

    # Combine the localizations
    smld_combined = combinelocalizations(smld_connected)

    return smld_connected, smld_combined
end
