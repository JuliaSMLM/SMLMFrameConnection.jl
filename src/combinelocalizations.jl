using SMLMData
using StatsBase

"""
    smld_combined = combinelocalizations(smld::BasicSMLD{T,E}) where {T, E<:SMLMData.AbstractEmitter}

Combine clustered localizations in `smld` into higher precision localizations.

# Description
This function combines localizations in `smld` that share the same value of
track_id.  Localizations are combined assuming they arose from
independent measurements of the same position with Gaussian errors.
"""
function combinelocalizations(smld::BasicSMLD{T,E}) where {T, E<:SMLMData.AbstractEmitter}
    # Extract arrays from emitters
    emitters = smld.emitters
    connectID = [e.track_id for e in emitters]
    sortindices = sortperm(connectID)
    connectID = connectID[sortindices]
    x = [e.x for e in emitters][sortindices]
    y = [e.y for e in emitters][sortindices]
    σ_x = [e.σ_x for e in emitters][sortindices]
    σ_y = [e.σ_y for e in emitters][sortindices]
    photons = [e.photons for e in emitters][sortindices]
    σ_photons = [e.σ_photons for e in emitters][sortindices]
    bg = [e.bg for e in emitters][sortindices]
    σ_bg = [e.σ_bg for e in emitters][sortindices]
    framenum = [e.frame for e in emitters][sortindices]
    datasetnum = [e.dataset for e in emitters][sortindices]
    id = [e.id for e in emitters][sortindices]

    # Loop over clusters and combine their localization data as appropriate.
    nperID = counts(connectID)
    nperID = nperID[nperID .!= 0]
    ncumulative = [0; cumsum(nperID)]
    nclusters = length(nperID)

    # Pre-allocate arrays for combined values
    # Use emitter's native precision, not SMLD's type parameter
    ET = typeof(first(emitters).x)
    x_combined = Vector{ET}(undef, nclusters)
    y_combined = Vector{ET}(undef, nclusters)
    σ_x_combined = Vector{ET}(undef, nclusters)
    σ_y_combined = Vector{ET}(undef, nclusters)
    photons_combined = Vector{ET}(undef, nclusters)
    σ_photons_combined = Vector{ET}(undef, nclusters)
    bg_combined = Vector{ET}(undef, nclusters)
    σ_bg_combined = Vector{ET}(undef, nclusters)
    connectID_combined = Vector{Int}(undef, nclusters)
    framenum_combined = Vector{Int}(undef, nclusters)
    datasetnum_combined = Vector{Int}(undef, nclusters)
    id_combined = Vector{Int}(undef, nclusters)

    for nn in 1:nclusters
        indices = (1:nperID[nn]) .+ ncumulative[nn]
        x_combined[nn] = StatsBase.mean(x[indices], weights(1.0 ./ σ_x[indices] .^ 2))
        y_combined[nn] = StatsBase.mean(y[indices], weights(1.0 ./ σ_y[indices] .^ 2))
        σ_x_combined[nn] = sqrt(1.0 / sum(1.0 ./ σ_x[indices] .^ 2))
        σ_y_combined[nn] = sqrt(1.0 / sum(1.0 ./ σ_y[indices] .^ 2))
        photons_combined[nn] = sum(photons[indices])
        σ_photons_combined[nn] = sqrt(sum(σ_photons[indices] .^ 2))
        bg_combined[nn] = sum(bg[indices])
        σ_bg_combined[nn] = sqrt(sum(σ_bg[indices] .^ 2))
        connectID_combined[nn] = connectID[indices[1]]
        framenum_combined[nn] = framenum[indices[1]]
        datasetnum_combined[nn] = datasetnum[indices[1]]
        id_combined[nn] = id[indices[1]]
    end

    # Create new emitters for combined results
    combined_emitters = Vector{SMLMData.Emitter2DFit{ET}}(undef, nclusters)
    for nn in 1:nclusters
        combined_emitters[nn] = SMLMData.Emitter2DFit{ET}(
            x_combined[nn], y_combined[nn],
            photons_combined[nn], bg_combined[nn],
            σ_x_combined[nn], σ_y_combined[nn],
            σ_photons_combined[nn], σ_bg_combined[nn],
            framenum_combined[nn], datasetnum_combined[nn],
            connectID_combined[nn], id_combined[nn]
        )
    end

    smld_combined = BasicSMLD(combined_emitters, smld.camera, smld.n_frames,
                              smld.n_datasets, copy(smld.metadata))

    return smld_combined
end
