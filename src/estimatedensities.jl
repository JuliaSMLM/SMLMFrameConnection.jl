using NearestNeighbors

"""
    estimatedensities(smld::BasicSMLD, clusterdata::Vector{Matrix{Float32}},
                      params::ParamStruct)

Estimate local emitter densities for clusters in `smld` and `clusterdata`.

# Description
The initial local densities `initialdensity` around each pre-cluster present in
`smld`/`clusterdata` are estimated based on the local density of pre-clusters
throughout the entire dataset as well as some of the rate parameters provided
in `params`.

# Returns
A vector of initial density estimates (emitters/pixel²) for each cluster.
"""
function estimatedensities(
    smld::BasicSMLD{T, E},
    clusterdata::Vector{Matrix{Float32}},
    params::ParamStruct
) where {T, E<:Emitter2DFit}
    # Calculate duty cycle
    dutycycle = params.k_on / (params.k_on + params.k_off + params.k_bleach)
    nclusters = length(clusterdata)

    # Determine maximum frame number
    maxframe = smld.n_frames > 0 ? smld.n_frames : maximum(e.frame for e in smld.emitters)

    # Handle single cluster case
    if nclusters == 1
        if size(clusterdata[1], 1) == 1
            return [1.0]
        else
            # Estimate density from cluster spatial extent
            x_data = clusterdata[1][:, 1]
            y_data = clusterdata[1][:, 2]
            clusterarea = (maximum(x_data) - minimum(x_data)) *
                         (maximum(y_data) - minimum(y_data))

            initialdensity = (1 / clusterarea) *
                            ((params.k_bleach / params.k_off) / (1 - params.p_miss)) /
                            (1 - exp(-params.k_bleach * dutycycle * (maxframe - 1)))
            return [initialdensity]
        end
    end

    # Compute cluster centers using inverse variance weighted means
    clustercenters = zeros(2, nclusters)
    for cluster_idx in 1:nclusters
        x = clusterdata[cluster_idx][:, 1]
        y = clusterdata[cluster_idx][:, 2]
        σ_x = clusterdata[cluster_idx][:, 3]
        σ_y = clusterdata[cluster_idx][:, 4]

        # MLE of position with Gaussian errors
        clustercenters[1, cluster_idx] = sum(x ./ σ_x.^2) / sum(1 ./ σ_x.^2)
        clustercenters[2, cluster_idx] = sum(y ./ σ_y.^2) / sum(1 ./ σ_y.^2)
    end

    # Estimate local cluster density from nearest neighbors
    kneighbors = min(params.nnearestclusters, nclusters - 1)
    kdtree = KDTree(clustercenters)
    _, nndist = knn(kdtree, clustercenters, kneighbors + 1)

    # Extract distances (excluding self at index 1)
    nndist = [nn[2:kneighbors+1] for nn in nndist]
    max_dist = [maximum(d) for d in nndist]

    # Cluster density: number of neighbors / circular area
    clusterdensity = (kneighbors + 1) ./ (π .* max_dist.^2)

    # Convert cluster density to underlying emitter density
    λ1 = params.k_bleach * dutycycle
    λ2 = (params.k_on + params.k_off + params.k_bleach) - λ1

    initialdensity = clusterdensity .*
                    (1.0 / dutycycle) * (1.0 / params.k_off) * (1.0 / (1.0 - params.p_miss)) ./
                    ((1.0 / λ1) * (1.0 - exp(-λ1 * (maxframe - 1.0))) -
                     (1.0 / λ2) * (1.0 - exp(-λ2 * (maxframe - 1.0))))

    return initialdensity
end
