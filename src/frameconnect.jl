using SMLMData

"""
    smld, smld_preclustered, smld_combined, params = frameconnect!(smld::SMLMData.SMLD2D; 
        nnearestclusters::Int=2, nsigmadev::Float64=5.0, maxframegap::Int=5, nmaxnn::Int=2)

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smld` are connected and
combined into higher precision localizations of that emitter.  This is done by
1) forming pre-clusters of localizations, 2) estimating rate parameters from
the pre-clusters, 3) solving a linear assignment problem for connecting 
localizations in each pre-cluster, and 4) combining the connected localizations
using their MLE position estimate assuming Gaussian noise.

# Inputs
-`smld`: SMLD2D structure containing the localizations that should be 
         connected (see SMLMData.SMLD2D for organization/fields).
-`nnearestclusters`: Number of nearest preclusters used for local density 
                     estimates. (default = 2)(see estimatedensities())
-`nsigmadev`: Multiplier of localization errors that defines a pre-clustering
              distance threshold. (default = 5)(see precluster())(pixels)
-`maxframegap`: Maximum frame gap between temporally adjacent localizations in
                a precluster. (default = 5)(see precluster())(frames)
-`nmaxnn`: Maximum number of nearest-neighbors inspected for precluster 
           membership.  Ideally, this would be set to inf, but that's not
           feasible for most data. (default = 2)(see precluster())

# Outputs
-`smld`: Input `smld` with field connectID updated to reflect connected
         localizations (however localizations remain uncombined).
-`smld_preclustered`: Copy of the input `smld` with the field connectID
                      populated to reflect the results of pre-clustering.
-`smld_combined`: Final frame-connection result (i.e., `smld` with 
                  localizations that seem to be from the same blinking event
                  combined into higher precision localizations).
-`params`: Structure of parameters used in the algorithm, with some copied 
           directly from the option kwargs to this function, and others 
           calculated internally (see SMLMFrameConnection.ParamStruct).
"""
function frameconnect!(smld::SMLMData.SMLD2D;
    nnearestclusters::Int = 2, nsigmadev::Float64 = 5.0,
    maxframegap::Int = 5, nmaxnn::Int = 2)

    # Prepare a ParamStruct to keep track of parameters used.
    params = ParamStruct()
    params.nnearestclusters = nnearestclusters
    params.nsigmadev = nsigmadev
    params.maxframegap = maxframegap
    params.nmaxnn = nmaxnn

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
    smld, smld_preclustered, smld_combined, params = frameconnect(smld::SMLMData.SMLD2D; 
        nnearestclusters::Int=2, nsigmadev::Float64=5.0, maxframegap::Int=5, nmaxnn::Int=2)

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smld` are connected and
combined into higher precision localizations of that emitter.  This is done by
1) forming pre-clusters of localizations, 2) estimating rate parameters from
the pre-clusters, 3) solving a linear assignment problem for connecting 
localizations in each pre-cluster, and 4) combining the connected localizations
using their MLE position estimate assuming Gaussian noise.

# Inputs
-`smld`: SMLD2D structure containing the localizations that should be 
         connected (see SMLMData.SMLD2D for organization/fields).
-`nnearestclusters`: Number of nearest preclusters used for local density 
                     estimates. (default = 2)(see estimatedensities())
-`nsigmadev`: Multiplier of localization errors that defines a pre-clustering
              distance threshold. (default = 5)(see precluster())(pixels)
-`maxframegap`: Maximum frame gap between temporally adjacent localizations in
                a precluster. (default = 5)(see precluster())(frames)
-`nmaxnn`: Maximum number of nearest-neighbors inspected for precluster 
           membership.  Ideally, this would be set to inf, but that's not
           feasible for most data. (default = 2)(see precluster())

# Outputs
-`smld`: Input `smld` with field connectID updated to reflect connected
         localizations (however localizations remain uncombined).
-`smld_preclustered`: Copy of the input `smld` with the field connectID
                      populated to reflect the results of pre-clustering.
-`smld_combined`: Final frame-connection result (i.e., `smld` with 
                  localizations that seem to be from the same blinking event
                  combined into higher precision localizations).
-`params`: Structure of parameters used in the algorithm, with some copied 
           directly from the option kwargs to this function, and others 
           calculated internally (see SMLMFrameConnection.ParamStruct).
"""
function frameconnect(smld::SMLMData.SMLD2D;
    nnearestclusters::Int = 2, nsigmadev::Float64 = 5,
    maxframegap::Int = 5, nmaxnn::Int = 2)

    return frameconnect!(deepcopy(smld);
        nnearestclusters = nnearestclusters, nsigmadev = nsigmadev,
        maxframegap = maxframegap, nmaxnn = nmaxnn)
end