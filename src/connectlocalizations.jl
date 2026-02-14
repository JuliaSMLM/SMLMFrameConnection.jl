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

    # Precompute unit-density rho_on curve (shared across all clusters).
    # Each cluster scales by its own rho_0.
    k_on = params.k_on
    k_off = params.k_off
    k_bleach = params.k_bleach
    dutycycle = k_on / (k_on + k_off + k_bleach)
    lambda1 = k_bleach * dutycycle
    lambda2 = (k_on + k_off + k_bleach) - lambda1
    base_rho_on = Vector{Float64}(undef, nframes)
    @inbounds for f in 1:nframes
        fm1 = Float64(f - 1)
        base_rho_on[f] = dutycycle * (exp(-lambda1 * fm1) - exp(-lambda2 * fm1))
    end

    # Loop through the clusters and solve the LAP.
    maxconnectID = maximum(connectID)
    idstocheck = findall(size.(clusterdata, 1) .> 1)
    for nn in idstocheck
        # Construct the cost matrix.
        costmatrix = create_costmatrix(clusterdata, params, nn, nframes, base_rho_on)

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