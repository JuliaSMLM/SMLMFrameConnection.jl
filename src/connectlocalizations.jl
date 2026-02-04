using StatsBase

"""
    connectlocalizations!(connectID::Vector{Int64},
        clusterdata::Vector{Matrix{Float32}},
        params::ParamStruct, nframes::Int64)

Connect localizations in `clusterdata` by solving a linear assignment problem.

# Description
Connect localizations in `clusterdata` by giving "connected" localizations the
same integer value for their field `smd.connectID`.  Associations are made by
solving a linear assignment problem which designates the ideal connections
between localizations.
"""
function connectlocalizations!(connectID::Vector{Int64},
    clusterdata::Vector{Matrix{Float32}}, params::ParamStruct, nframes::Int64)

    # Loop through the clusters and solve the LAP.
    maxconnectID = maximum(connectID)
    idstocheck = findall(size.(clusterdata, 1) .> 1)
    for nn in idstocheck
        # Construct the cost matrix.
        costmatrix = create_costmatrix(clusterdata, params, nn, nframes)

        # Solve the LAP.
        assignment = solveLAP(costmatrix)

        # Update `connectID` to reflect the LAP solution.
        _, maxconnectID = linkclusters!(connectID, maxconnectID,
            Int64.(clusterdata[nn][:, 9]), assignment[1])  # sortindex is column 9
    end

    # Ensure `connectID` is contains the set of integers 1:nclusters.
    compress_connectID!(connectID)

    return connectID
end

"""
    connectID = connectlocalizations(connectID::Vector{Int64}, 
        clusterdata::Vector{Matrix{Float32}}, 
        params::ParamStruct, nframes::Int64)

Connect localizations in `clusterdata` by solving a linear assignment problem.

# Description
Connect localizations in `clusterdata` by giving "connected" localizations the
same integer value for their field `smd.connectID`.  Associations are made by
solving a linear assignment problem which designates the ideal connections
between localizations.
"""
function connectlocalizations(connectID::Vector{Int64},
    clusterdata::Vector{Matrix{Float32}}, params::ParamStruct, nframes::Int64)

    return connectlocalizations!(deepcopy(connectID), clusterdata, params, nframes)
end