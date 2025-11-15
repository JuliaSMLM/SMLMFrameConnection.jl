"""
    frameconnect!(smld::BasicSMLD; nnearestclusters::Int=2, nsigmadev::Float64=5.0,
                  maxframegap::Int=5, nmaxnn::Int=2)

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smld` are connected and
combined into higher precision localizations. This is done via a multi-step process:

1. Form pre-clusters of localizations based on spatial and temporal proximity
2. Estimate rate parameters (k_on, k_off, k_bleach, p_miss) from the pre-clusters
3. Solve a linear assignment problem (LAP) for connecting localizations within each pre-cluster
4. Combine connected localizations using MLE position estimates assuming Gaussian noise

# Arguments
- `smld`: BasicSMLD containing the localizations to be connected
- `nnearestclusters`: Number of nearest preclusters for local density estimates (default = 2)
- `nsigmadev`: Multiplier of localization errors for pre-clustering distance threshold (default = 5)
- `maxframegap`: Maximum frame gap for temporally adjacent localizations in a precluster (default = 5)
- `nmaxnn`: Maximum number of nearest-neighbors for precluster membership (default = 2)

# Returns
- `smld`: Input SMLD with `track_id` field updated to reflect connections
- `smld_preclustered`: SMLD with `track_id` values from pre-clustering stage
- `smld_combined`: Final result with connected localizations combined into higher precision localizations
- `params`: ParamStruct containing all algorithm parameters (including estimated ones)
"""
function frameconnect!(
    smld::BasicSMLD{T, E};
    nnearestclusters::Int = 2,
    nsigmadev::Float64 = 5.0,
    maxframegap::Int = 5,
    nmaxnn::Int = 2
) where {T, E<:Emitter2DFit}
    # Generate pre-clusters of localizations
    smld_preclustered = precluster(smld, ParamStruct(
        Float64[],
        nnearestclusters,
        0.0, 0.0, 0.0, 0.0,
        nsigmadev,
        maxframegap,
        nmaxnn
    ))

    # Organize clusters for analysis
    clusterdata = organizeclusters(smld_preclustered)

    # Estimate rate parameters from pre-clusters
    k_on, k_off, k_bleach, p_miss, _ = estimateparams(smld_preclustered, clusterdata)

    # Create params struct with estimated values
    params_temp = ParamStruct(
        Float64[],
        nnearestclusters,
        k_on, k_off, k_bleach, p_miss,
        nsigmadev,
        maxframegap,
        nmaxnn
    )

    # Estimate underlying density of emitters
    initialdensity = estimatedensities(smld_preclustered, clusterdata, params_temp)

    # Create final params struct
    params = ParamStruct(
        initialdensity,
        nnearestclusters,
        k_on, k_off, k_bleach, p_miss,
        nsigmadev,
        maxframegap,
        nmaxnn
    )

    # Connect localizations by solving the LAP
    n_frames = smld.n_frames > 0 ? smld.n_frames : maximum(e.frame for e in smld.emitters)

    # Get track_ids from preclustering and update via LAP
    precluster_track_ids = [e.track_id for e in smld_preclustered.emitters]
    connected_track_ids = connectlocalizations(precluster_track_ids, clusterdata, params, n_frames)

    # Create new SMLD with updated track_ids
    new_emitters = [set_track_id(smld.emitters[i], connected_track_ids[i]) for i in 1:length(smld.emitters)]
    smld_connected = BasicSMLD(new_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)

    # Combine connected localizations into higher precision localizations
    smld_combined = combinelocalizations(smld_connected)

    return smld_connected, smld_preclustered, smld_combined, params
end

"""
    frameconnect(smld::BasicSMLD; nnearestclusters::Int=2, nsigmadev::Float64=5.0,
                 maxframegap::Int=5, nmaxnn::Int=2)

Non-mutating version of `frameconnect!`.

# Description
Connects repeated localizations of the same emitter without modifying the input.
See `frameconnect!` for detailed documentation.

# Arguments
- `smld`: BasicSMLD containing the localizations to be connected
- `nnearestclusters`: Number of nearest preclusters for local density estimates (default = 2)
- `nsigmadev`: Multiplier of localization errors for pre-clustering distance threshold (default = 5)
- `maxframegap`: Maximum frame gap for temporally adjacent localizations in a precluster (default = 5)
- `nmaxnn`: Maximum number of nearest-neighbors for precluster membership (default = 2)

# Returns
- `smld`: New SMLD with `track_id` field updated to reflect connections
- `smld_preclustered`: New SMLD with `track_id` values from pre-clustering stage
- `smld_combined`: Final result with connected localizations combined
- `params`: ParamStruct containing all algorithm parameters
"""
function frameconnect(
    smld::BasicSMLD{T, E};
    nnearestclusters::Int = 2,
    nsigmadev::Float64 = 5.0,
    maxframegap::Int = 5,
    nmaxnn::Int = 2
) where {T, E<:Emitter2DFit}
    return frameconnect!(
        deepcopy(smld);
        nnearestclusters=nnearestclusters,
        nsigmadev=nsigmadev,
        maxframegap=maxframegap,
        nmaxnn=nmaxnn
    )
end
