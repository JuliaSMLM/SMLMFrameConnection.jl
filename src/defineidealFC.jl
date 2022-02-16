using SMLMData


"""
    connect1DS!(smld::SMLMData.SMLD2D, dataset::Int; maxframegap::Int = 5)

Define the "ideal" frame-connection result for simulated `smld` with one dataset.

# Description
This is a helper function which defines the ideal FC result for a single 
dataset.  The user should probably not bother using this, and instead should
call defineidealFC().

# Inputs
- `smld`: SMLMData.SMLD2D with smld.connectID populated to indicate emitter ID.
- `dataset`: Dataset number to be connected.
- `maxframegap`: Maximum frame gap allowed between localizations connected in
                 the "ideal" result.
"""
function connect1DS!(smld::SMLMData.SMLD2D, dataset::Int; maxframegap::Int = 5)
    # Loop through associated localizations and combine them as appropriate.
    currentDS = smld.datasetnum .== dataset
    smld_current = smld[currentDS]
    maxID = maximum(smld.connectID)
    smld_current.connectID .+= maxID # recompressed later
    for ii in unique(smld_current.connectID)
        # Determine which localizations belong to the ii-th emitter.
        currentinds = findall(smld_current.connectID .== ii)
    
        # Sort these localizations with respect to their frame.
        framenum = smld_current.framenum[currentinds]
        sortinds = sortperm(framenum)
        sortedframes = framenum[sortinds]
        sorted_currentinds = currentinds[sortinds]
    
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
                smld_current.connectID[sorted_currentinds[ff+1]] =
                    smld_current.connectID[sorted_currentinds[ff]]
            else
                # Localization ff+1 should be a new blinking event.
                maxID += 1
                smld_current.connectID[sorted_currentinds[ff+1]] = maxID
            end
        end
    end

    # Update smld.connectID for the current dataset.
    smld.connectID[currentDS] = smld_current.connectID

    return smld
end

function connect1DS(smld::SMLMData.SMLD2D, dataset::Int; maxframegap::Int = 5)
    return connect1DS!(deepcopy(smld))
end

"""
    smld, smld_combined = defineidealFC!(smld::SMLMData.SMLD2D; 
        maxframegap::Int = 5)

Define the "ideal" frame-connection result for a simulated `smld`.

# Description
This function defines the "ideal" frame connection result from a simulation.
That is to say, for a simulated SMLD2D structure `smld` with `connectID` field
populated to indicate emitter membership of localizations, this function will
generate an "ideal" FC result which combines all blinking events that appeared
with frame gaps less than `maxframegap` of one another.  Note that for very
high duty cycles, multiple blinking events might be mistakingly combined by
this method (i.e., if the emitter blinks back on within `maxframegap` frames
of its previous blink).  Note that localizations are not allowed to be 
connected across datasets.

# Inputs
-`smld`: SMLMData.SMLD2D with smld.connectID populated to indicate emitter ID.
-`maxframegap`: Maximum frame gap allowed between localizations connected in
                the "ideal" result.

# Outputs
-`smld_out`: Ideal frame-connection result on input `smld` with localizations
             combined as appropriate.
-`smld_connected`: Copy of the input `smld` with smld_connected.connectID 
                   modified to reflect blinking event ID (i.e., `smld_out` is
                   generated as smld_out = combinelocalizations(smld_connected))
"""
function defineidealFC!(smld::SMLMData.SMLD2D; maxframegap::Int = 5)
    # Loop through datasets and combine localizations as appropriate.
    for ii in unique(smld.datasetnum)
        connect1DS!(smld, ii; maxframegap = maxframegap)
    end
    compress_connectID!(smld.connectID)

    return smld, combinelocalizations(smld)
end

function defineidealFC(smld::SMLMData.SMLD2D; maxframegap::Int = 5)
    return defineidealFC!(deepcopy(smld); maxframegap = maxframegap)
end