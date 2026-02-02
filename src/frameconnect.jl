using SMLMData

"""
    (combined, info) = frameconnect(smld::BasicSMLD{T,E}; kwargs...)

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smld` are connected and
combined into higher precision localizations of that emitter. This is done by
1) forming pre-clusters of localizations, 2) estimating rate parameters from
the pre-clusters, 3) solving a linear assignment problem for connecting
localizations in each pre-cluster, and 4) combining the connected localizations
using their MLE position estimate assuming Gaussian noise.

# Arguments
- `smld::BasicSMLD`: Localizations to connect. Must contain emitters with valid
                     position uncertainties (σ_x, σ_y).

# Keyword Arguments
- `nnearestclusters::Int=2`: Number of nearest preclusters used for local density
                             estimates (see `estimatedensities`)
- `nsigmadev::Float64=5.0`: Multiplier of localization errors that defines a
                            pre-clustering distance threshold (see `precluster`)
- `maxframegap::Int=5`: Maximum frame gap between temporally adjacent localizations
                        in a precluster (see `precluster`)
- `nmaxnn::Int=2`: Maximum number of nearest-neighbors inspected for precluster
                   membership (see `precluster`)

# Returns
A tuple `(combined, info)`:
- `combined::BasicSMLD`: Connected localizations combined into higher precision results
- `info::ConnectInfo`: Track assignments and algorithm metadata (see [`ConnectInfo`](@ref))

# Example
```julia
(combined, info) = frameconnect(smld)
println("Connected \$(info.n_input) → \$(info.n_combined) localizations")
println("Formed \$(info.n_tracks) tracks from \$(info.n_preclusters) preclusters")

# Access track assignments for downstream analysis
track_ids = [e.track_id for e in info.connected.emitters]
```
"""
function frameconnect(smld::BasicSMLD{T,E};
    nnearestclusters::Int = 2, nsigmadev::Float64 = 5.0,
    maxframegap::Int = 5, nmaxnn::Int = 2) where {T, E<:SMLMData.AbstractEmitter}

    t_start = time_ns()

    # Prepare a ParamStruct to keep track of parameters used.
    params = ParamStruct()
    params.nnearestclusters = nnearestclusters
    params.nsigmadev = nsigmadev
    params.maxframegap = maxframegap
    params.nmaxnn = nmaxnn

    # Generate pre-clusters of localizations in `smld`.
    smld_preclustered = precluster(smld, params)
    clusterdata = organizeclusters(smld_preclustered)
    n_preclusters = length(clusterdata)

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

    elapsed_ns = time_ns() - t_start

    # Build ConnectInfo
    n_tracks = length(unique(connectID_final))
    info = ConnectInfo{T}(
        smld_connected,
        length(smld.emitters),
        n_tracks,
        length(smld_combined.emitters),
        params.k_on,
        params.k_off,
        params.k_bleach,
        params.p_miss,
        params.initialdensity,
        elapsed_ns,
        :lap,
        n_preclusters
    )

    return (smld_combined, info)
end
