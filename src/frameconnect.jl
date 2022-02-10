using SMLMData

"""
    smld, smld_preclustered, smld_combined, params = frameconnect!(smld::SMLMData.SMLD2D, 
        params::ParamStruct)

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smd` are connected and
combined into higher precision localizations of that emitter.  This is done by
1) forming pre-clusters of localizations, 2) estimating rate parameters from
the pre-clusters, 3) solving a linear assignment problem for connecting 
localizations in each pre-cluster, and 4) combining the connected localizations
using their MLE position estimate assuming Gaussian noise.
"""
function frameconnect!(smld::SMLMData.SMLD2D, params::ParamStruct = ParamStruct())
    # Make a copy of `params` (some values are updated within this function,
    # but the user may need to keep their original structure).
    params = deepcopy(params)

    # Generate pre-clusters of localizations in `smld`.
    smld_preclustered = precluster(smld, params)
    clusterdata = organizeclusters(smld_preclustered)

    # Estimate rate parameters.
    params.k_on, params.k_off, params.k_bleach, params.p_miss =
        estimateparams(smld_preclustered, clusterdata)

    # Estimate the underlying density of emitters.
    params.initialdensity =
        estimatedensities(smld_preclustered, clusterdata, params)

    # Connect localizations in `smld` by solving the LAP.
    nframes = isempty(smld.nframes) ? maximum(smld.framenum) : smld.nframes
    smld.connectID = connectlocalizations(smld_preclustered.connectID,
        clusterdata, params, nframes)

    # Combine the connected localizations into higher precision localizations.
    smld_combined = combinelocalizations(smld)

    return smld, smld_preclustered, smld_combined, params
end

"""
    smld, smld_preclustered, smld_combined, params = frameconnect(smld::SMLMData.SMLD2D, 
        params::ParamStruct)

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smd` are connected and
combined into higher precision localizations of that emitter.  This is done by
1) forming pre-clusters of localizations, 2) estimating rate parameters from
the pre-clusters, 3) solving a linear assignment problem for connecting 
localizations in each pre-cluster, and 4) combining the connected localizations
using their MLE position estimate assuming Gaussian noise.
"""
function frameconnect(smld::SMLMData.SMLD2D, params::ParamStruct = ParamStruct())
    return frameconnect!(deepcopy(smld), params)
end