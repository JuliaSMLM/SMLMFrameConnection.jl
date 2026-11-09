using SMLMData
using NearestNeighbors
using Statistics

"""
    initial_density = estimatedensities(smld::BasicSMLD{T,E},
        clusterdata::Vector{Matrix{Float64}}, params::ParamStruct) where {T, E<:SMLMData.AbstractEmitter}

Estimate local emitter densities for clusters in `smld` and `clusterdata`.

# Description
The initial local densities `initial_density` around each pre-cluster present in
`smld`/`clusterdata` are estimated based on the local density of pre-clusters
throughout the entire set of data as well as some of the rate parameters
provided in `params`.
"""
function estimatedensities(smld::BasicSMLD{T,E},
    clusterdata::Vector{Matrix{Float64}},
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
            # Localizations sharing an x or a y value span zero area; floor the area at
            # the localization error ellipse so the density stays finite.
            locarea = pi * Statistics.mean(clusterdata[1][:, 3]) *
                      Statistics.mean(clusterdata[1][:, 4])
            clusterarea = max(clusterarea, locarea)
            initial_density = (1 / clusterarea) *
                             ((params.k_bleach / params.k_off) / (1 - params.p_miss)) /
                             (1 - exp(-params.k_bleach * dutycycle * (maxframe - 1)))
        end
        return [initial_density]  # Return as vector to match ParamStruct.initial_density type
    end

    # Determine the center of all clusters assuming each arose from the same
    # emitter using precision-weighted mean with full covariance.
    clustercenters = zeros(2, nclusters)
    clusterscale = zeros(nclusters)   # localization error scale √(mean σ_x · mean σ_y)
    for nn = 1:nclusters
        # Column layout: 1:x 2:y 3:σ_x 4:σ_y 5:σ_xy 6:frame ...
        # Precision-weighted mean with full covariance: center = (Σ Pᵢ)⁻¹ Σ Pᵢ posᵢ
        c = clusterdata[nn]
        μx, μy, _, _, _ = _precision_mean(view(c, :, 1), view(c, :, 2), view(c, :, 3),
                                          view(c, :, 4), view(c, :, 5), axes(c, 1))
        clustercenters[1, nn] = μx
        clustercenters[2, nn] = μy
        clusterscale[nn] = sqrt(Statistics.mean(view(c, :, 3)) * Statistics.mean(view(c, :, 4)))
    end

    # Estimate the local cluster density based on the distance to
    # nearest-neighbors.
    kneighbors = minimum([params.n_density_neighbors; nclusters - 1])
    kneighbors = max(kneighbors, 1)  # Ensure at least 1 neighbor
    # Ensure we don't request more neighbors than available (including self)
    kneighbors = min(kneighbors, nclusters - 1)

    kdtree = NearestNeighbors.KDTree(clustercenters)
    _, nndist_raw = NearestNeighbors.knn(kdtree, clustercenters, kneighbors + 1)

    # Skip the cluster itself and any cluster at the same center, then take the k-th
    # nearest remaining neighbor. "The same center" means closer than a millionth of the
    # localization error: no measurement resolves that, and a center a rounding error
    # away from an equal one (a different summation order gives one) would otherwise
    # count as a neighbor at 1e-15 μm and a density of 1e29 per μm².
    nndist = zeros(nclusters)
    for i in 1:nclusters
        same = 1e-6 * clusterscale[i]
        nonzero_dists = filter(d -> d > same, nndist_raw[i])
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
