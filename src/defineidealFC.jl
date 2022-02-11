using SMLMData


"""
# Description
This function defines the "ideal" frame connection result from a simulation.
That is to say, for a simulated SMLD2D structure `smld` with `connectID` field
populated to indicate emitter membership of localizations, this function will
generate an "ideal" FC result which combines all blinking events that appeared
with frame gaps less than `maxframegap` of one another.  Note that for very
high duty cycles, multiple blinking events might be mistakingly combined by
this method (i.e., if the emitter blinks back on within `maxframegap` frames
of its previous blink).
"""
function defineidealFC(smld::SMLMData.SMLD2D; maxframegap::Int = 5)
    # Loop through associated localizations and combine them as appropriate.
    smld_out = SMLMData.SMLD2D(0)
    smld_connected = deepcopy(smld)
    maxID = maximum(smld.connectID)
    for ii in unique(smld.connectID)
        # Determine which localizations belong to the ii-th emitter.
        currentinds = findall(smld.connectID .== ii)

        # Sort these localizations with respect to their frame.
        framenum = smld.framenum[currentinds]
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
                smld_connected.connectID[sorted_currentinds[ff + 1]] =
                    smld_connected.connectID[sorted_currentinds[ff]]
            else
                # Localization ff+1 should be a new blinking event.
                maxID += 1
                smld_connected.connectID[sorted_currentinds[ff+1]] = maxID
            end
        end
    end
    compress_connectID!(smld_connected.connectID)
    smld_out = combinelocalizations(smld_connected)

    return smld_out, smld_connected
end



# """
# # Description
# This function defines the "ideal" frame connection result from a simulation.
# That is to say, for a simulated SMLD2D structure `smld` with `connectID` field
# populated to indicate emitter membership of localizations, this function will
# generate an "ideal" FC result which combines all blinking events that appeared
# with frame gaps less than `maxframegap` of one another.  Note that for very
# high duty cycles, multiple blinking events might be mistakingly combined by
# this method (i.e., if the emitter blinks back on within `maxframegap` frames
# of its previous blink).
# """
# function defineidealFC(smld::SMLMData.SMLD2D; maxframegap::Int = 5)
#     # Loop through associated localizations and combine them as appropriate.
#     smld_out = SMLMData.SMLD2D(0)
#     for ii in unique(smld.connectID)
#         # Determine which localizations belong to the ii-th emitter.
#         currentinds = findall(smld.connectID .== ii)

#         # Sort these localizations with respect to their frame.
#         framenum = smld.framenum[currentinds]
#         sortinds = sortperm(framenum)
#         sortedframes = framenum[sortinds]
#         sorted_currentinds = currentinds[sortinds]

#         # Check if all localizations are within the framegap, combining them
#         # immediately if that is the case.
#         # println(diff(sortedframes))
#         if all(diff(sortedframes) .<= maxframegap)
#             smld_connected = deepcopy(smld)
#             smld_combined = combinelocalizations(smld[currentinds])
#             smld_out = SMLMData.catsmld(smld_out, smld_combined)
#             continue
#         end

#         # Determine which localizations we can combine.
#         smld_connected = deepcopy(smld[sorted_currentinds])
#         for ff = 1:(length(sortedframes)-1)
#             if (sortedframes[ff+1]-sortedframes[ff]) <= maxframegap
#                 # Connect these localizations.
#                 smld_connected.connectID[[ff; ff+1]] .= 
#                     minimum(smld_connected.connectID[[ff; ff+1]])
#             end
#         end
#         smld_combined = combinelocalizations(smld_connected)
#         smld_connectedALL = SMLMData.catsmld(smld_connectedALL, smld_connected)
#         smld_out = SMLMData.catsmld(smld_out, smld_combined)
#     end

#     return smld_out, smld_connected
# end