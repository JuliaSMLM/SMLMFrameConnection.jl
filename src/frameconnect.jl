using SMLMData

"""
    (combined, info) = frameconnect(smld::BasicSMLD, config::FrameConnectConfig)
    (combined, info) = frameconnect(smld::BasicSMLD; kwargs...)

Connect repeated localizations of the same emitter in `smld`.

# Description
Repeated localizations of the same emitter present in `smld` are connected and
combined into higher precision localizations of that emitter. This is done by
1) forming pre-clusters of localizations, 2) estimating rate parameters from
the pre-clusters, 3) solving a linear assignment problem for connecting
localizations in each pre-cluster, and 4) combining the connected localizations
using their MLE position estimate assuming Gaussian noise.

# Arguments
- `smld::BasicSMLD`: 2D localizations to connect, in any order. Must contain emitters
                     with valid position uncertainties (σ_x, σ_y). Localizations that
                     carry z (3D) are refused with an `ArgumentError`; output emitters are
                     `Emitter2DFit`, so fields beyond those (e.g. a fit p-value) are not
                     carried into `info.connected`.
- `config::FrameConnectConfig`: Configuration parameters (optional, can use kwargs instead)

# Keyword Arguments (equivalent to FrameConnectConfig fields)
- `n_density_neighbors::Int=2`: Number of nearest preclusters used for local density
                             estimates (see `estimatedensities`)
- `max_sigma_dist::Float64=5.0`: Multiplier of localization errors that defines a
                            pre-clustering distance threshold (see `precluster`)
- `max_frame_gap::Int=5`: Maximum frame gap between temporally adjacent localizations
                        in a precluster (see `precluster`)
- `max_neighbors::Int=2`: Maximum number of nearest-neighbors inspected for precluster
                   membership (see `precluster`)
- `track_length::Union{Tuple{Float64,Float64},Nothing}=nothing`: Inclusive `(min, max)`
                   range on localizations per track; tracks outside it are dropped. `nothing`
                   is a no-op. `(2.0, Inf)` removes single-frame blinks; `(2.0, 50.0)` also
                   drops over-long tracks (fiducials, sticky strands). The number dropped is
                   reported as `info.n_filtered` (see [`combinelocalizations`](@ref)).

# Returns
A tuple `(combined, info)`:
- `combined::BasicSMLD`: Connected localizations combined into higher precision results
- `info::FrameConnectInfo`: Track assignments and algorithm metadata (see [`FrameConnectInfo`](@ref))

# Example
```julia
# Using kwargs (most common)
(combined, info) = frameconnect(smld)
(combined, info) = frameconnect(smld; max_frame_gap=10)

# Using config struct
config = FrameConnectConfig(max_frame_gap=10, max_sigma_dist=3.0)
(combined, info) = frameconnect(smld, config)

println("Connected \$(info.n_input) → \$(info.n_combined) localizations")
println("Formed \$(info.n_tracks) tracks from \$(info.n_preclusters) preclusters")

# Access track assignments for downstream analysis
track_ids = [e.track_id for e in info.connected.emitters]
```
"""
function frameconnect(smld::BasicSMLD{T,E}; kwargs...) where {T, E<:SMLMData.AbstractEmitter}
    # kwargs form forwards to config form
    config = FrameConnectConfig(; kwargs...)
    return frameconnect(smld, config)
end

function frameconnect(smld::BasicSMLD{T,E}, config::FrameConnectConfig) where {T, E<:SMLMData.AbstractEmitter}
    t_start = time()
    _check_2d(smld, "frameconnect")

    # Degrade gracefully on empty input: there is nothing to connect, so return
    # the (empty) input unchanged with a zeroed info. This short-circuits before
    # the downstream pipeline (precluster, estimateparams, ...), which otherwise
    # reduces over empty collections and throws (e.g. `counts` on an empty frame
    # vector, or `zero(Type{Any})` for an untyped empty emitter vector).
    if isempty(smld.emitters)
        # Mirror the non-empty path's calibration contract: when calibration is
        # configured, still emit a (graceful fallback) CalibrationResult rather
        # than `nothing`, so downstream code that reads `info.calibration.*` sees
        # the same shape for 0 and 1-2 localizations. analyze_calibration is
        # empty-safe (0 pairs -> calibration_applied=false fallback).
        cal_result = config.calibration === nothing ? nothing :
                     analyze_calibration(smld, config.calibration)
        info = FrameConnectInfo{T}(
            smld,                 # connected: input passthrough (no tracks to assign)
            0, 0, 0, 0,           # n_input, n_tracks, n_combined, n_filtered
            0.0, 0.0, 0.0, 0.0,   # k_on, k_off, k_bleach, p_miss
            Float64[],            # initial_density
            time() - t_start,     # elapsed_s
            :lap,                 # algorithm
            0,                    # n_preclusters
            cal_result,           # calibration: nothing if disabled, else fallback
        )
        return (smld, info)
    end

    # Run the pipeline on a copy sorted by (dataset, frame, x, y), so the result does
    # not depend on the order of the input; track IDs are mapped back to input order.
    perm = sortperm(smld.emitters; by = e -> (e.dataset, e.frame, e.x, e.y))
    smld_sorted = BasicSMLD(smld.emitters[perm], smld.camera, smld.n_frames,
                            smld.n_datasets, smld.metadata)

    # Prepare a ParamStruct to keep track of parameters used.
    params = ParamStruct()
    params.n_density_neighbors = config.n_density_neighbors
    params.max_sigma_dist = config.max_sigma_dist
    params.max_frame_gap = config.max_frame_gap
    params.max_neighbors = config.max_neighbors

    # Generate pre-clusters of localizations in `smld`.
    smld_preclustered = precluster(smld_sorted, params)
    clusterdata = organizeclusters(smld_preclustered)
    n_preclusters = length(clusterdata)
    _warn_large_preclusters(clusterdata, params)

    # Estimate rate parameters.
    params.k_on, params.k_off, params.k_bleach, params.p_miss =
        estimateparams(smld_preclustered, clusterdata)

    # Estimate the underlying density of emitters.
    params.initial_density =
        estimatedensities(smld_preclustered, clusterdata, params)

    # Get nframes
    nframes = smld.n_frames > 0 ? smld.n_frames : maximum(e.frame for e in smld.emitters)

    # Connect localizations in `smld` by solving the LAP.
    # Extract track_id from preclustered emitters
    connectID_precluster = [e.track_id for e in smld_preclustered.emitters]
    connectID_sorted = connectlocalizations(connectID_precluster,
        clusterdata, params, nframes)
    connectID_final = Vector{Int}(undef, length(connectID_sorted))
    connectID_final[perm] = connectID_sorted

    # Connected localizations, in the sorted order: calibration and combination run on
    # these, so every stage sees one order and `combined` does not depend on the input
    # order either. `info.connected` is the same localizations in the caller's order.
    # Use emitter's native precision, not SMLD's type parameter
    emitters_sorted = smld_sorted.emitters
    ET = typeof(first(emitters_sorted).x)  # Get precision from emitter fields
    connected_sorted = Vector{SMLMData.Emitter2DFit{ET}}(undef, length(emitters_sorted))
    for i in eachindex(emitters_sorted)
        e = emitters_sorted[i]
        connected_sorted[i] = SMLMData.Emitter2DFit{ET}(
            e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
            e.frame, e.dataset, connectID_sorted[i], e.id
        )
    end
    smld_connected_sorted = BasicSMLD(connected_sorted, smld.camera, smld.n_frames,
                                      smld.n_datasets, copy(smld.metadata))
    connected_caller = similar(connected_sorted)
    connected_caller[perm] = connected_sorted
    smld_connected = BasicSMLD(connected_caller, smld.camera, smld.n_frames,
                               smld.n_datasets, copy(smld.metadata))

    # Optional uncertainty calibration (between connection and combination)
    cal_result = nothing
    smld_to_combine = smld_connected_sorted
    if config.calibration !== nothing
        cal_result = analyze_calibration(smld_connected_sorted, config.calibration)
        if cal_result.calibration_applied
            smld_to_combine = apply_calibration(smld_connected_sorted, cal_result)
        end
    end

    # Combine the connected localizations into higher precision localizations.
    smld_combined = combinelocalizations(smld_to_combine;
        track_length=config.track_length)

    elapsed_s = time() - t_start

    # Build FrameConnectInfo
    n_tracks = length(unique(connectID_final))

    # Tracks dropped specifically by the length filter (cause-attributed count,
    # not n_tracks - n_combined, so it stays correct if combine ever drops
    # clusters for other reasons). connectID_final is compressed to 1:n_tracks.
    n_filtered = if config.track_length === nothing
        0
    else
        lo, hi = config.track_length
        count(n -> !(lo <= n <= hi), counts(connectID_final))
    end
    info = FrameConnectInfo{T}(
        smld_connected,
        length(smld.emitters),
        n_tracks,
        length(smld_combined.emitters),
        n_filtered,
        params.k_on,
        params.k_off,
        params.k_bleach,
        params.p_miss,
        params.initial_density,
        elapsed_s,
        :lap,
        n_preclusters,
        cal_result
    )

    return (smld_combined, info)
end

# The LAP on a precluster of N localizations is a 2N x 2N assignment problem, whose time
# grows as N^3 (about a minute at N = 2000). Preclusters are the connected components of
# the neighbor links, so very dense data (several active emitters per um^2 per frame) or
# dense clusters can join into very large ones.
const _LARGE_PRECLUSTER = 2000

function _warn_large_preclusters(clusterdata, params::ParamStruct)
    isempty(clusterdata) && return nothing
    sizes = size.(clusterdata, 1)
    nmax = maximum(sizes)
    if nmax > _LARGE_PRECLUSTER
        @warn "frameconnect: $(count(>(_LARGE_PRECLUSTER), sizes)) precluster(s) hold more " *
              "than $(_LARGE_PRECLUSTER) localizations (largest: $nmax). Linking time grows " *
              "as the cube of a precluster's size; a smaller max_sigma_dist " *
              "($(params.max_sigma_dist)) or max_frame_gap ($(params.max_frame_gap)) " *
              "keeps preclusters smaller."
    end
    return nothing
end
