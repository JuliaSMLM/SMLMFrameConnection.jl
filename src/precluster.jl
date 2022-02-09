using SMLMData
using Statistics
using NearestNeighbors

"""
    smld_preclustered = precluster(smld::SMLMData.SMLD2D, 
        params::ParamStruct = ParamStruct())

Cluster localizations in `smld` based on distance and time thresholds in `params`.

# Description
Localizations in the input structure `smd` are clustered together based on
their spatiotemporal separations.  All localizations within a spatial 
threshold of `params.nsigmadev*mean([smld.σ_x smld.σ_y)` and a temporal 
threshold of `params.maxframegap` of one another will be clustered together,
meaning that these localizations now share the same unique integer value for
smld.connectID.

# Notes
Pre-clustering allows localizations observed in the same frame to be in the
same cluster.  This is done to prevent exclusion of the "correct" localization
from its ideal cluster due to a previously included "incorrect" localization
into that cluster.
"""
function precluster(smld::SMLMData.SMLD2D, params::ParamStruct = ParamStruct())
    # Make a copy of smld, grab some of its fields (to improve speed), and
    # sort w.r.t. framenum.
    sortindicesDS = sortperm(smld.datasetnum)
    framenum = smld.framenum[sortindicesDS]
    datasetnum = smld.datasetnum[sortindicesDS]
    xy = transpose([smld.x smld.y])
    xy = xy[:, sortindicesDS]
    σ_x = smld.σ_x[sortindicesDS]
    σ_y = smld.σ_y[sortindicesDS]
    mean_se = Statistics.mean([σ_x σ_y]; dims = 2)

    # Isolate some parameters from params.
    maxframegap = params.maxframegap
    nsigmadev = params.nsigmadev

    # Initialize a connectID array, with each localization being considered
    # a unique cluster.
    nlocalizations = length(framenum)
    connectID = collect(1:nlocalizations)

    # Loop through frames and add localizations to clusters.
    nperdataset = counts(datasetnum)
    nperdataset = nperdataset[nperdataset.!=0]
    ncumulativeDS = [0; cumsum(nperdataset)]
    datasets = unique(datasetnum)
    maxID = nlocalizations
    for nn = 1:length(datasets)
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
        frames = unique(framenumCDS)
        clusterinds = [[ind] for ind in 1:nperdataset[nn]]
        for ff = 1:length(frames)
            # Determine which localizations should be considered for
            # clustering.
            currentindFN = (1:nperframe[ff]) .+ ncumulativeFN[ff]
            candidateind = findall((framenumCDS .>= (frames[ff].-maxframegap)) .& 
                                   (framenumCDS .<= frames[ff]))
            if length(candidateind) < 2
                maxID += 1
                connectIDCDS[currentindFN] = (1:nperframe[ff]) .+ maxID
                maxID += nperframe[ff]
                continue
            end

            # Find the nearest-neighbor to the current localizations
            # which is a candidate for clustering.
            kdtree = KDTree(xyCDS[:, candidateind])
            nnindices, nndist = knn(kdtree, xyCDS[:, currentindFN],
                min(params.nmaxnn + 1, length(candidateind)), true)

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

                for jj in updateinds
                    clusterinds[jj] = [clusterinds[jj]; updateinds]
                end
            end
        end
        connectID[currentindDS] = connectIDCDS
    end

    # Store the updated connectID in the output structure, ensuring that the
    # values range from 1:NClusters (which is expected by later codes).
    smld_preclustered = deepcopy(smld)
    smld_preclustered.connectID = compress_connectID(connectID[sortindicesDS])

    return smld_preclustered
end