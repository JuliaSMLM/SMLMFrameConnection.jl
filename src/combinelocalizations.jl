using SMLMData
using StatsBase

"""
    smld_combined = combinelocalizations(smld::SMLMData.SMLD2D)

Combine clustered localizations in `smld` into higher precision localizations.

# Description
This function combines localizations in `smld` that share the same value of
smld.connectID.  Localizations are combined assuming they arose from 
independent measurements of the same position with Gaussian errors.
"""
function combinelocalizations(smld::SMLMData.SMLD2D)
    # Isolate and sort some arrays from `smld`.
    connectID = smld.connectID
    sortindices = sortperm(connectID)
    connectID = connectID[sortindices]
    x = smld.x[sortindices]
    y = smld.y[sortindices]
    σ_x = smld.σ_x[sortindices]
    σ_y = smld.σ_y[sortindices]
    photons = smld.photons[sortindices]
    bg = smld.bg[sortindices]
    framenum = smld.framenum[sortindices]
    datasetnum = smld.datasetnum[sortindices]

    # Loop over clusters and combine their localization data as appropriate.
    nperID = counts(connectID)
    nperID = nperID[nperID .!= 0]
    ncumulative = [0; cumsum(nperID)]
    nclusters = length(nperID)
    smld_combined = deepcopy(smld)
    loopinds = (1:nclusters)[nperID .> 1]
    for nn in loopinds
        indices = (1:nperID[nn]) .+ ncumulative[nn]
        smld_combined.x[nn] = StatsBase.mean(x[indices], weights(1.0 / σ_x[indices] .^ 2))
        smld_combined.y[nn] = StatsBase.mean(y[indices], weights(1.0 / σ_y[indices] .^ 2))
        smld_combined.σ_x[nn] = sqrt(1.0 / sum(1.0 ./ σ_x[indices] .^ 2))
        smld_combined.σ_y[nn] = sqrt(1.0 / sum(1.0 ./ σ_y[indices] .^ 2))
        smld_combined.photons[nn] = sum(photons[indices])
        smld_combined.bg[nn] = sum(bg[indices])
        smld_combined.connectID[nn] = connectID[indices[1]]
        smld_combined.framenum[nn] = framenum[indices[1]]
        smld_combined.datasetnum[nn] = datasetnum[indices[1]]
    end
    smld_combined = SMLMData.isolatesmld(smld_combined, 1:nclusters)

    return smld_combined
end