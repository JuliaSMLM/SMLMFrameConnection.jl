"""
    connectID, maxconnectID = linkclusters!(connectID::Int64, 
        maxconnectID::Int64, updateind::Vector{Int64}, 
        assignment::Vector{Int64})

Force localizations linked by `assignment` to share same the same `connectID`
"""
function linkclusters!(connectID::Vector{Int64}, maxconnectID::Int64,
    updateind::Union{Vector{Int64},Int64}, assignment::Vector{Int64})

    # Initialize each localization being updated as a new cluster.
    nlocalizations = div(length(assignment), 2)
    connectID[updateind] = maxconnectID .+ collect(1:nlocalizations)

    # Loop through the first `nlocalizations` entries of assignment and make
    # the appropriate assignments in `connectID`.
    # NOTE: The cost matrix organization means that the last `nlocalizations`
    #       entries of assignment don't require changes in `connectID`.
    connectID_current = connectID[updateind]
    connectID_copy = deepcopy(connectID_current)
    for nn = 1:nlocalizations
        # If this is a cluster "birth", no action is needed.
        if assignment[nn] > nlocalizations
            continue
        end

        # Determine which ID should be attributed with this assignment.
        newID = minimum([connectID_current[nn]
            connectID_current[assignment[nn]]
            connectID_copy[nn]
            connectID_copy[assignment[nn]]])

        # Update the reassigned values of `connectID`.
        connectID_current[[nn assignment[nn]]] .= newID
    end

    # Store the updated `connectID`s ensuring that they are all greater than
    # the `maxconnectID`.
    connectID[updateind] = connectID_current .+ maxconnectID
    maxconnectID += nlocalizations

    return connectID, maxconnectID
end

"""
    connectID, maxconnectID = linkclusters(
        connectID::Int64, maxconnectID::Int64, 
        updateind::Vector{Int64}, assignment::Vector{Int64})

Force localizations linked by `assignment` to share same the same `connectID`
"""
function linkclusters(connectID::Vector{Int64}, maxconnectID::Int64,
    updateind::Union{Vector{Int64},Int64}, assignment::Vector{Int64})

    return linkclusters!(deepcopy(connectID), deepcopy(maxconnectID), updateind, assignment)
end