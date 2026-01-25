# Test fixtures and helper functions for SMLMFrameConnection tests

using SMLMData

"""
    make_test_camera(; xpixels=512, ypixels=512, pixelsize=0.1)

Create an IdealCamera for testing.
"""
function make_test_camera(; xpixels=512, ypixels=512, pixelsize=0.1)
    return IdealCamera(1:xpixels, 1:ypixels, pixelsize)
end

"""
    make_emitter(x, y, frame; σ_xy=0.02, photons=1000.0, track_id=0)

Create a single Emitter2DFit for testing with sensible defaults.
"""
function make_emitter(x::Real, y::Real, frame::Int;
                      σ_xy::Real=0.02,
                      photons::Real=1000.0,
                      bg::Real=10.0,
                      track_id::Int=0,
                      dataset::Int=1)
    return Emitter2DFit{Float64}(
        Float64(x), Float64(y),
        Float64(photons), Float64(bg),
        Float64(σ_xy), Float64(σ_xy),
        Float64(sqrt(photons)), Float64(sqrt(bg)),
        frame, dataset, track_id, 0
    )
end

"""
    make_test_smld(emitters; n_frames=nothing, n_datasets=1)

Create a BasicSMLD from a vector of emitters.
"""
function make_test_smld(emitters::Vector{<:Emitter2DFit};
                        n_frames::Union{Int,Nothing}=nothing,
                        n_datasets::Int=1)
    cam = make_test_camera()
    nf = isnothing(n_frames) ? maximum(e.frame for e in emitters) : n_frames
    return BasicSMLD(emitters, cam, nf, n_datasets)
end

"""
    make_blinking_molecule(x, y, frames; σ_xy=0.02, position_jitter=0.002)

Create emitters representing a single molecule blinking across multiple frames.
Uses deterministic small offsets to simulate localization uncertainty.
"""
function make_blinking_molecule(x::Real, y::Real, frames::AbstractVector{Int};
                                σ_xy::Real=0.02,
                                position_jitter::Real=0.002,
                                track_id::Int=0)
    emitters = Emitter2DFit{Float64}[]
    for (i, f) in enumerate(frames)
        # Deterministic small offset pattern (not random for reproducibility)
        offset = position_jitter * ((i % 3) - 1)  # -jitter, 0, +jitter pattern
        push!(emitters, make_emitter(x + offset, y + offset, f; σ_xy=σ_xy, track_id=track_id))
    end
    return emitters
end

"""
    make_two_molecules_smld(; separation=1.0)

Create test SMLD with two well-separated molecules, each blinking for 3 frames.
"""
function make_two_molecules_smld(; separation::Real=1.0)
    emitters = vcat(
        make_blinking_molecule(5.0, 5.0, [1, 2, 3]; track_id=1),
        make_blinking_molecule(5.0 + separation, 5.0, [1, 2, 3]; track_id=2)
    )
    return make_test_smld(emitters; n_frames=3)
end

"""
    make_single_molecule_smld(; n_frames=5)

Create test SMLD with a single molecule blinking across n_frames consecutive frames.
"""
function make_single_molecule_smld(; n_frames::Int=5, x::Real=5.0, y::Real=5.0)
    emitters = make_blinking_molecule(x, y, collect(1:n_frames); track_id=1)
    return make_test_smld(emitters; n_frames=n_frames)
end

"""
    make_prelabeled_smld()

Create SMLD with track_id pre-populated (for testing combinelocalizations).
"""
function make_prelabeled_smld()
    emitters = [
        # Cluster 1: two localizations
        make_emitter(5.0, 5.0, 1; track_id=1),
        make_emitter(5.01, 5.01, 2; track_id=1),
        # Cluster 2: three localizations
        make_emitter(10.0, 10.0, 1; track_id=2),
        make_emitter(10.01, 9.99, 2; track_id=2),
        make_emitter(9.99, 10.01, 3; track_id=2),
        # Cluster 3: single localization
        make_emitter(15.0, 15.0, 1; track_id=3),
    ]
    return make_test_smld(emitters; n_frames=3)
end
