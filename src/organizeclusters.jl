using SMLMData

"""
    clusterdata = organizeclusters(smld::SMLMData.SMLD2D)

Organize pre-clusters into a vector indexing distinct clusters.

# Description
Pre-clusters in the input `smld`, as related by their shared integer value k_off
`smld.connectID`, are organized into a vector of matrices `clusterdata`.  Each
index of `clusterdata` corresponds to a distinct cluster.
"""
function organizeclusters(smld::SMLMData.SMLD2D)
    # Isolate and sort some arrays from the smd structure.
    connectID = smld.connectID
    sortindices = sortperm(connectID)
    combineddata = [smld.x smld.y smld.σ_x smld.σ_y smld.framenum smld.datasetnum smld.connectID]
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