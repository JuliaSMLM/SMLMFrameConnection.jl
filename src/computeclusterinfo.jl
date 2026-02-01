"""
    clusterdurations, nobservations = 
        computeclusterinfo(clusterdata::Vector{Matrix{Float32}})

Compute the durations of and number of member localizations in clusters.

# Description
This method computes the duration of each cluster and the number of
localizations in each cluster present in the input `clusterdata`.
"""
function computeclusterinfo(clusterdata::Vector{Matrix{Float32}})
    # Compute some quantities needed from each of the preclusters.
    nclusters = size(clusterdata)[1]
    clusterdurations = ones(nclusters)
    nobservations = ones(nclusters)
    for nn = 1:nclusters
        # Compute the duration of this cluster.
        # Column layout: 1:x 2:y 3:σ_x 4:σ_y 5:σ_xy 6:frame 7:dataset 8:connectID 9:sortindex
        currentframes = clusterdata[nn][:, 6]
        clusterdurations[nn] = maximum(currentframes) - minimum(currentframes) + 1

        # Compute the number of observations of this cluster excluding the 
        # multiple observations of localizations in the same frame (that is,
        # ignore the generous pre-clustering artifact of multiple localizations
        # per frame, since this calculation assumes we have just one emitter).
        nobservations[nn] = length(unique(currentframes))
    end

    return clusterdurations, nobservations
end