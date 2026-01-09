using SMLMData

"""
    clusterdata = organizeclusters(smld::BasicSMLD{T,SMLMData.Emitter2DFit{T}}) where T

Organize pre-clusters into a vector indexing distinct clusters.

# Description
Pre-clusters in the input `smld`, as related by their shared integer value of
`track_id`, are organized into a vector of matrices `clusterdata`.  Each
index of `clusterdata` corresponds to a distinct cluster.
"""
function organizeclusters(smld::BasicSMLD{T,SMLMData.Emitter2DFit{T}}) where T
    # Extract arrays from emitters
    emitters = smld.emitters
    connectID = [e.track_id for e in emitters]
    x = [e.x for e in emitters]
    y = [e.y for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    framenum = [e.frame for e in emitters]
    datasetnum = [e.dataset for e in emitters]

    # Sort by connectID
    sortindices = sortperm(connectID)
    combineddata = [x y σ_x σ_y framenum datasetnum connectID]
    combineddata = combineddata[sortindices, :]

    # Loop over preclusters and organize member localizations into a more
    # accessible format.
    npercluster = counts(connectID)
    npercluster = npercluster[npercluster.!=0]
    ncumulative = [0; cumsum(npercluster)]
    uniqueIDs = unique(connectID)
    nclusters = length(uniqueIDs)
    clusterdata = Vector{Matrix{Float32}}(undef, nclusters)
    for nn = 1:nclusters
        clusterind = (1:npercluster[nn]) .+ ncumulative[nn]
        clusterdata[nn] = [combineddata[clusterind, :] sortindices[clusterind]]
    end

    return clusterdata
end
