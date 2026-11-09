using SMLMData
using Statistics
using NearestNeighbors

"""
    smld_preclustered = precluster(smld::BasicSMLD{T,E},
        params::ParamStruct = ParamStruct()) where {T, E<:SMLMData.AbstractEmitter}

Cluster localizations in `smld` based on distance and time thresholds in `params`.

# Description
Localizations in the input structure `smld` are clustered together based on
their spatiotemporal separations. Each localization is linked to those of its
`params.max_neighbors` nearest neighbors, among the localizations of the same
dataset in its own frame and the `params.max_frame_gap` frames before it, that lie
within `params.max_sigma_dist * (s_i + s_j)`, where `s` is a localization's mean of
σ_x and σ_y. Preclusters are the connected components of these links: their
localizations share the same unique integer value for their track_id field. The
result does not depend on the order of the input.

# Notes
Pre-clustering allows localizations observed in the same frame to be in the
same cluster.  This is done to prevent exclusion of the "correct" localization
from its ideal cluster due to a previously included "incorrect" localization
into that cluster.
"""
function precluster(smld::BasicSMLD{T,E},
                    params::ParamStruct = ParamStruct()) where {T, E<:SMLMData.AbstractEmitter}
    # Extract arrays from emitters
    emitters = smld.emitters
    n_emitters = length(emitters)

    # Nothing to cluster: return the (empty) input unchanged. The reductions
    # below (`counts` on datasetnum/framenum, `first(emitters)`) assume at least
    # one localization and otherwise throw on an empty collection.
    n_emitters == 0 && return smld

    # Work in (dataset, frame) order (stable), so the clusters do not depend on the
    # order of the input.
    framenum = [e.frame for e in emitters]
    datasetnum = [e.dataset for e in emitters]
    order = sortperm(collect(zip(datasetnum, framenum)))
    fr = framenum[order]
    ds = datasetnum[order]
    ET = typeof(first(emitters).x)
    xy = Matrix{ET}(undef, 2, n_emitters)
    mean_se = Vector{ET}(undef, n_emitters)
    for (k, i) in enumerate(order)
        e = emitters[i]
        xy[1, k] = e.x
        xy[2, k] = e.y
        mean_se[k] = (e.σ_x + e.σ_y) / 2
    end

    # Isolate some parameters from params.
    max_frame_gap = params.max_frame_gap
    max_sigma_dist = params.max_sigma_dist
    max_neighbors = params.max_neighbors

    # Link each localization to those of its `max_neighbors` nearest neighbors, among
    # the localizations of its own frame and the `max_frame_gap` frames before it, that
    # lie within the distance threshold. Clusters are the connected components of
    # these links (union-find; each root is the smallest index in its component).
    parent = collect(1:n_emitters)
    lo = 1
    while lo <= n_emitters
        hi = lo
        while hi < n_emitters && ds[hi+1] == ds[lo]
            hi += 1
        end
        frames_ds = view(fr, lo:hi)
        a = lo
        while a <= hi
            b = a
            while b < hi && fr[b+1] == fr[a]
                b += 1
            end
            c = lo - 1 + searchsortedfirst(frames_ds, fr[a] - max_frame_gap)
            ncand = b - c + 1
            if ncand >= 2
                kdtree = NearestNeighbors.KDTree(xy[:, c:b])
                nnindices, nndist = NearestNeighbors.knn(kdtree, xy[:, a:b],
                    min(max_neighbors + 1, ncand), true)
                for (ii, i) in enumerate(a:b)
                    for (jj, d) in zip(nnindices[ii], nndist[ii])
                        j = c - 1 + jj
                        if d <= max_sigma_dist * (mean_se[i] + mean_se[j])
                            _uf_union!(parent, i, j)
                        end
                    end
                end
            end
            a = b + 1
        end
        lo = hi + 1
    end

    # Label by component root, compress to 1:NClusters, and return to input order.
    connectID = [_uf_find!(parent, k) for k in 1:n_emitters]
    connectID_compressed = Vector{Int}(undef, n_emitters)
    connectID_compressed[order] = compress_connectID(connectID)

    # Create new emitters with track_id set
    # Use emitter's native precision, not SMLD's type parameter
    new_emitters = Vector{SMLMData.Emitter2DFit{ET}}(undef, n_emitters)
    for i in 1:n_emitters
        e = emitters[i]
        new_emitters[i] = SMLMData.Emitter2DFit{ET}(
            e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
            e.frame, e.dataset, connectID_compressed[i], e.id
        )
    end

    smld_preclustered = BasicSMLD(new_emitters, smld.camera, smld.n_frames,
                                   smld.n_datasets, copy(smld.metadata))

    return smld_preclustered
end

# Union-find with path halving; the root of a component is its smallest index.
function _uf_find!(parent::Vector{Int}, i::Int)
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    return i
end

function _uf_union!(parent::Vector{Int}, i::Int, j::Int)
    ri = _uf_find!(parent, i)
    rj = _uf_find!(parent, j)
    ri == rj && return
    if ri < rj
        parent[rj] = ri
    else
        parent[ri] = rj
    end
    return
end
