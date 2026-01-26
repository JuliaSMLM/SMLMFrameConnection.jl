using SMLMData

"""
    result = frameconnect(smld::BasicSMLD{T,E};
        nnearestclusters::Int=2, nsigmadev::Float64=5.0,
        maxframegap::Int=5, nmaxnn::Int=2) where {T, E<:SMLMData.AbstractEmitter}

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smld` are connected and
combined into higher precision localizations of that emitter.  This is done by
1) forming pre-clusters of localizations, 2) estimating rate parameters from
the pre-clusters, 3) solving a linear assignment problem for connecting
localizations in each pre-cluster, and 4) combining the connected localizations
using their MLE position estimate assuming Gaussian noise.

# Inputs
- `smld`: BasicSMLD containing the localizations that should be connected.
          Must contain emitters with valid position uncertainties (σ_x, σ_y).
- `nnearestclusters`: Number of nearest preclusters used for local density
                      estimates. (default = 2)(see estimatedensities())
- `nsigmadev`: Multiplier of localization errors that defines a pre-clustering
               distance threshold. (default = 5)(see precluster())(microns)
- `maxframegap`: Maximum frame gap between temporally adjacent localizations in
                 a precluster. (default = 5)(see precluster())(frames)
- `nmaxnn`: Maximum number of nearest-neighbors inspected for precluster
            membership.  Ideally, this would be set to inf, but that's not
            feasible for most data. (default = 2)(see precluster())

# Outputs
Returns a NamedTuple with fields:
- `combined`: Final frame-connection result (i.e., `smld` with localizations that
              seem to be from the same blinking event combined into higher
              precision localizations).
- `connected`: Input `smld` with track_id field updated to reflect connected
               localizations (localizations remain uncombined).
- `params`: Structure of parameters used in the algorithm, with some copied
            directly from the option kwargs to this function, and others
            calculated internally (see SMLMFrameConnection.ParamStruct).
"""
function frameconnect(smld::BasicSMLD{T,E};
    nnearestclusters::Int = 2, nsigmadev::Float64 = 5.0,
    maxframegap::Int = 5, nmaxnn::Int = 2) where {T, E<:SMLMData.AbstractEmitter}

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

    # Get nframes
    nframes = smld.n_frames > 0 ? smld.n_frames : maximum(e.frame for e in smld.emitters)

    # Connect localizations in `smld` by solving the LAP.
    # Extract track_id from preclustered emitters
    connectID_precluster = [e.track_id for e in smld_preclustered.emitters]
    connectID_final = connectlocalizations(connectID_precluster,
        clusterdata, params, nframes)

    # Create smld_connected with updated track_id
    # Use emitter's native precision, not SMLD's type parameter
    emitters = smld.emitters
    ET = typeof(first(emitters).x)  # Get precision from emitter fields
    new_emitters = Vector{SMLMData.Emitter2DFit{ET}}(undef, length(emitters))
    for i in 1:length(emitters)
        e = emitters[i]
        new_emitters[i] = SMLMData.Emitter2DFit{ET}(
            e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
            e.frame, e.dataset, connectID_final[i], e.id
        )
    end
    smld_connected = BasicSMLD(new_emitters, smld.camera, smld.n_frames,
                                smld.n_datasets, copy(smld.metadata))

    # Combine the connected localizations into higher precision localizations.
    smld_combined = combinelocalizations(smld_connected)

    return (combined=smld_combined, connected=smld_connected, params=params)
end
