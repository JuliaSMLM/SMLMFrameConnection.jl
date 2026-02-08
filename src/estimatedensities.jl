using SMLMData
using NearestNeighbors

"""
    initial_density = estimatedensities(smld::BasicSMLD{T,E},
        clusterdata::Vector{Matrix{Float32}}, params::ParamStruct) where {T, E<:SMLMData.AbstractEmitter}

Estimate local emitter densities for clusters in `smld` and `clusterdata`.

# Description
The initial local densities `initial_density` around each pre-cluster present in
`smld`/`clusterdata` are estimated based on the local density of pre-clusters
throughout the entire set of data as well as some of the rate parameters
provided in `params`.
"""
function estimatedensities(smld::BasicSMLD{T,E},
    clusterdata::Vector{Matrix{Float32}},
    params::ParamStruct) where {T, E<:SMLMData.AbstractEmitter}

    # Define some new parameters.
    dutycycle = params.k_on / (params.k_on + params.k_off + params.k_bleach)
    nclusters = length(clusterdata)

    # Get max frame from n_frames or from emitters
    if smld.n_frames > 0
        maxframe = smld.n_frames
    else
        maxframe = maximum(e.frame for e in smld.emitters)
    end

    # If only one cluster is present, we should return an answer right away and stop.
    if nclusters == 1
        if size(clusterdata[1], 1) == 1
            initial_density = 1.0
        else
            clusterarea = (maximum(clusterdata[1][:, 1]) - minimum(clusterdata[1][:, 1])) *
                          (maximum(clusterdata[1][:, 2]) - minimum(clusterdata[1][:, 2]))
            initial_density = (1 / clusterarea) *
                             ((params.k_bleach / params.k_off) / (1 - params.p_miss)) /
                             (1 - exp(-params.k_bleach * dutycycle * (maxframe - 1)))
        end
        return [initial_density]  # Return as vector to match ParamStruct.initial_density type
    end

    # Determine the center of all clusters assuming each arose from the same
    # emitter using precision-weighted mean with full covariance.
    clustercenters = zeros(2, nclusters)
    for nn = 1:nclusters
        # Isolate some arrays to improve readability.
        # Column layout: 1:x 2:y 3:σ_x 4:σ_y 5:σ_xy 6:frame ...
        x = clusterdata[nn][:, 1]
        y = clusterdata[nn][:, 2]
        σ_x = clusterdata[nn][:, 3]
        σ_y = clusterdata[nn][:, 4]
        σ_xy = clusterdata[nn][:, 5]

        # Compute the cluster centers using precision-weighted mean with full covariance.
        # P = Σ⁻¹ = [σ_y² -σ_xy; -σ_xy σ_x²] / det where det = σ_x²σ_y² - σ_xy²
        # Combined: center = (Σ Pᵢ)⁻¹ * Σ(Pᵢ * posᵢ)
        P_sum = zeros(2, 2)
        μ_weighted = zeros(2)
        for i in eachindex(x)
            det_i = σ_x[i]^2 * σ_y[i]^2 - σ_xy[i]^2
            det_i = max(det_i, eps())  # Ensure positive definite
            # Precision matrix elements
            P11 = σ_y[i]^2 / det_i
            P22 = σ_x[i]^2 / det_i
            P12 = -σ_xy[i] / det_i
            P_sum[1,1] += P11
            P_sum[2,2] += P22
            P_sum[1,2] += P12
            P_sum[2,1] += P12
            μ_weighted[1] += P11 * x[i] + P12 * y[i]
            μ_weighted[2] += P12 * x[i] + P22 * y[i]
        end
        # Solve P_sum * center = μ_weighted
        clustercenters[1:2, nn] = P_sum \ μ_weighted
    end

    # Estimate the local cluster density based on the distance to
    # nearest-neighbors.
    kneighbors = minimum([params.n_density_neighbors; nclusters - 1])
    kneighbors = max(kneighbors, 1)  # Ensure at least 1 neighbor
    # Ensure we don't request more neighbors than available (including self)
    kneighbors = min(kneighbors, nclusters - 1)

    kdtree = NearestNeighbors.KDTree(clustercenters)
    _, nndist_raw = NearestNeighbors.knn(kdtree, clustercenters, kneighbors + 1)

    # Get non-zero distances (skip self-distance which is 0)
    # Then get the k-th nearest non-self neighbor
    nndist = zeros(nclusters)
    for i in 1:nclusters
        # Filter out zero distances (self)
        nonzero_dists = filter(d -> d > 0, nndist_raw[i])
        if length(nonzero_dists) >= kneighbors
            sort!(nonzero_dists)
            nndist[i] = nonzero_dists[kneighbors]  # k-th nearest non-self neighbor
        elseif !isempty(nonzero_dists)
            nndist[i] = maximum(nonzero_dists)
        else
            # Fallback: use a reasonable distance in microns to avoid Inf
            nndist[i] = 1.0
        end
    end

    clusterdensity = (kneighbors + 1) ./ (pi * nndist .^ 2)

    # Estimate the density of underlying emitters based on cluster density.
    lambda1 = params.k_bleach * dutycycle
    lambda2 = (params.k_on + params.k_off + params.k_bleach) - lambda1
    initial_density = clusterdensity *
                     (1.0 / dutycycle) * (1.0 / params.k_off) * (1.0 / (1.0 - params.p_miss)) ./
                     ((1.0 / lambda1) * (1.0 - exp(-lambda1 * (maxframe - 1.0))) -
                      (1.0 / lambda2) * (1.0 - exp(-lambda2 * (maxframe - 1.0))))

    return initial_density
end
