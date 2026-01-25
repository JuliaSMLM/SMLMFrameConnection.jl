using SMLMData
using Statistics
using NearestNeighbors

"""
    smld_preclustered = precluster(smld::BasicSMLD{T,E},
        params::ParamStruct = ParamStruct()) where {T, E<:SMLMData.AbstractEmitter}

Cluster localizations in `smld` based on distance and time thresholds in `params`.

# Description
Localizations in the input structure `smld` are clustered together based on
their spatiotemporal separations.  All localizations within a spatial
threshold of `params.nsigmadev*mean([σ_x σ_y])` and a temporal
threshold of `params.maxframegap` of one another will be clustered together,
meaning that these localizations now share the same unique integer value for
their track_id field.

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

    framenum = [e.frame for e in emitters]
    datasetnum = [e.dataset for e in emitters]
    x = [e.x for e in emitters]
    y = [e.y for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]

    # Sort w.r.t. datasetnum
    sortindicesDS = sortperm(datasetnum)
    framenum = framenum[sortindicesDS]
    datasetnum = datasetnum[sortindicesDS]
    xy = transpose([x y])
    xy = xy[:, sortindicesDS]
    σ_x = σ_x[sortindicesDS]
    σ_y = σ_y[sortindicesDS]
    mean_se = Statistics.mean([σ_x σ_y]; dims = 2)

    # Isolate some parameters from params.
    maxframegap = params.maxframegap
    nsigmadev = params.nsigmadev
    nmaxnn = params.nmaxnn

    # Initialize a connectID array, with each localization being considered
    # a unique cluster.
    nlocalizations = length(framenum)
    connectID = collect(1:nlocalizations)

    # Loop through frames and add localizations to clusters.
    nperdataset = StatsBase.counts(datasetnum)
    nperdataset = nperdataset[nperdataset.!=0]
    ncumulativeDS = [0; cumsum(nperdataset)]
    maxID = nlocalizations
    for nn = 1:length(unique(datasetnum))
        # Isolate some arrays for the nn-th dataset and sort w.r.t. their
        # framenum.
        currentindDS = (1:nperdataset[nn]) .+ ncumulativeDS[nn]
        sortindicesFN = sortperm(framenum[currentindDS])
        currentindDS = currentindDS[sortindicesFN]
        framenumCDS = framenum[currentindDS]
        xyCDS = xy[:, currentindDS]
        mean_seCDS = mean_se[currentindDS]
        connectIDCDS = connectID[currentindDS]

        # Loop through frames and add localizations to clusters.
        nperframe = counts(framenumCDS)
        nperframe = nperframe[nperframe.!=0]
        ncumulativeFN = [0; cumsum(nperframe)]
        clusterinds = [[ind] for ind in 1:nperdataset[nn]]
        frames = unique(framenumCDS)
        for ff = 1:length(frames)
            # Determine which localizations should be considered for
            # clustering.
            currentindFN = (1:nperframe[ff]) .+ ncumulativeFN[ff]
            candidateind = findall((framenumCDS .>= (frames[ff] .- maxframegap)) .&
                                   (framenumCDS .<= frames[ff]))
            if length(candidateind) < 2
                maxID += 1
                connectIDCDS[currentindFN] = (1:nperframe[ff]) .+ maxID
                maxID += nperframe[ff]
                continue
            end

            # Find the nearest-neighbor to the current localizations
            # which is a candidate for clustering.
            kdtree = NearestNeighbors.KDTree(xyCDS[:, candidateind])
            nnindices, nndist = NearestNeighbors.knn(kdtree, xyCDS[:, currentindFN],
                min(nmaxnn + 1, length(candidateind)), true)

            # Assign localizations to clusters based on `nndist`.
            for ii in 1:nperframe[ff]
                # Determine which candidates meet our distance cutoff.
                se_sum = mean_seCDS[currentindFN[ii]] .+
                         mean_seCDS[candidateind[nnindices[ii]]]
                validnninds = nnindices[ii][nndist[ii].<=(nsigmadev*se_sum)]

                # Update connectIDCDS to reflect the new clusters.
                updateinds = unique([currentindFN[ii]
                    candidateind[validnninds]
                    findall(connectIDCDS .== connectIDCDS[currentindFN[ii]])
                    clusterinds[candidateind[validnninds]][1]])
                connectIDCDS[updateinds] .= minimum(connectIDCDS[updateinds])

                # Keep track of which indices have been associated.
                for jj in updateinds
                    clusterinds[jj] = [clusterinds[jj][1]; updateinds]
                end
            end
        end
        connectID[currentindDS] = deepcopy(connectIDCDS)
    end

    # Compress connectID to range from 1:NClusters
    connectID_compressed = Vector{Int}(undef, nlocalizations)
    connectID_compressed[sortindicesDS] = compress_connectID(connectID)

    # Create new emitters with track_id set
    # Use emitter's native precision, not SMLD's type parameter
    ET = typeof(first(emitters).x)
    new_emitters = Vector{SMLMData.Emitter2DFit{ET}}(undef, n_emitters)
    for i in 1:n_emitters
        e = emitters[i]
        new_emitters[i] = SMLMData.Emitter2DFit{ET}(
            e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg,
            e.frame, e.dataset, connectID_compressed[i], e.id
        )
    end

    smld_preclustered = BasicSMLD(new_emitters, smld.camera, smld.n_frames,
                                   smld.n_datasets, copy(smld.metadata))

    return smld_preclustered
end
