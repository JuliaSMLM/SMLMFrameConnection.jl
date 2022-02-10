using StatsBase

"""
    compress_connectID!(connectID::Vector{Int32})

Make `connectID` consists of the set of integers `1:length(unique(connectID))``
"""
function compress_connectID!(connectID::Vector{Int64})
    # Force `connectID` to consist of the integers 
    # `1:length(unique(connectID))` without skipping any integers.
    sortindices = sortperm(connectID)
    connectID_sorted = deepcopy(connectID[sortindices])
    connectIDlist = unique(connectID_sorted)
    nperID = counts(connectID_sorted)
    nperID = nperID[nperID.!=0]
    ncumulative = [0; cumsum(nperID)]
    for nn = 1:length(connectIDlist)
        clusterind = (1:nperID[nn]) .+ ncumulative[nn]
        connectID_sorted[clusterind] .= nn
    end
    connectID[sortindices] = deepcopy(connectID_sorted)

    return connectID
end

"""
    connectID = compress_connectID(connectID::Vector{Int32})

Make `connectID` consists of the set of integers `1:length(unique(connectID))``
"""
function compress_connectID(connectID::Vector{Int64})
    return compress_connectID!(deepcopy(connectID))
end