using SMLMData
using NearestNeighbors

"""
    initialdensity = estimatedensities(smld::SMLMData.SMLD2D, 
        clusterdata::Vector{Matrix{Float32}}, params::ParamStruct)

Estimate local emitter densities for clusters in `smld` and `clusterdata`.

# Description
The initial local densities `initialdensity` around each pre-cluster present in
`smld`/`clusterdata` are estimated based on the local density of pre-clusters
throughout the entire set of data as well as some of the rate parameters 
provided in `params`.
"""
function estimatedensities(smld::SMLMData.SMLD2D,
    clusterdata::Vector{Matrix{Float32}},
    params::ParamStruct)

    # Define some new parameters.
    dutycycle = params.k_on / (params.k_on + params.k_off + params.k_bleach)
    nclusters = length(clusterdata)
    if isempty(smld.nframes)
        maxframe = maximum(smld.framenum)
    else
        maxframe = smld.nframes
    end

    # If only one cluster is present, we should return an answer right away and stop.
    if nclusters == 1
        if size(clusterdata[1], 1) == 1
            initialdensity = 1
        else
            clusterarea = (maximum(clusterdata[1][:, 1]) - minimum(clusterdata[1][:, 1])) *
                          (maximum(clusterdata[1][:, 1]) - minimum(clusterdata[1][:, 1]))
            initialdensity = (1 / clusterarea) *
                             ((params.k_bleach / params.k_off) / (1 - params.p_miss)) /
                             (1 - exp(-params.k_bleach * dutycycle * (maxframe - 1)))
        end
        return initialdensity
    end

    # Determine the center of all clusters assuming each arose from the same
    # emitter.
    clustercenters = zeros(2, nclusters)
    for nn = 1:nclusters
        # Isolate some arrays to improve readability.
        x = clusterdata[nn][:, 1]
        y = clusterdata[nn][:, 2]
        σ_x = clusterdata[nn][:, 3]
        σ_y = clusterdata[nn][:, 4]

        # Compute the cluster centers based on MLE of position.
        clustercenters[1:2, nn] = [sum(x ./ σ_x .^ 2) / sum(1 ./ σ_x .^ 2)
            sum(y ./ σ_y .^ 2) / sum(1 ./ σ_y .^ 2)]
    end

    # Estimate the local cluster density based on the distance to 
    # nearest-neighbors.
    kneighbors = minimum([params.nnearestclusters nclusters - 1])
    kdtree = NearestNeighbors.KDTree(clustercenters)
    _, nndist = NearestNeighbors.knn(kdtree, clustercenters, kneighbors + 1)
    nndist = getindex.(nndist, 1:(kneighbors-1))
    clusterdensity = (kneighbors + 1) ./ (pi * nndist .^ 2)

    # Estimate the density of underlying emitters based on cluster density.
    lambda1 = params.k_bleach * dutycycle
    lambda2 = (params.k_on + params.k_off + params.k_bleach) - lambda1
    initialdensity = clusterdensity *
                     (1.0 / dutycycle) * (1.0 / params.k_off) * (1.0 / (1.0 - params.p_miss)) ./
                     ((1.0 / lambda1) * (1.0 - exp(-lambda1 * (maxframe - 1.0))) -
                      (1.0 / lambda2) * (1.0 - exp(-lambda2 * (maxframe - 1.0))))

    return initialdensity
end