# Example script for basic frame connection usage.
#
# This example demonstrates:
# 1. Creating synthetic SMLM data with blinking molecules
# 2. Running frame connection to combine repeated localizations
# 3. Comparing results with ground truth

using SMLMFrameConnection
using SMLMData
using Random
using Statistics

# Set random seed for reproducibility
Random.seed!(42)

#=
Generate synthetic blinking molecule data.

We simulate multiple molecules, each blinking across several frames with:
- Gaussian-distributed position noise (localization uncertainty)
- Variable photon counts
- Random blinking patterns
=#

function generate_synthetic_smld(;
    n_molecules::Int = 20,
    n_frames::Int = 100,
    fov_size::Float64 = 10.0,    # μm
    σ_loc::Float64 = 0.02,       # localization uncertainty (μm)
    p_on::Float64 = 0.3,         # probability of being "on" each frame
    mean_photons::Float64 = 1000.0
)
    emitters = Emitter2DFit{Float64}[]

    # Generate random molecule positions
    mol_x = fov_size * rand(n_molecules)
    mol_y = fov_size * rand(n_molecules)

    for mol_id in 1:n_molecules
        true_x, true_y = mol_x[mol_id], mol_y[mol_id]

        for frame in 1:n_frames
            # Random blinking - molecule is "on" with probability p_on
            if rand() < p_on
                # Add localization noise
                obs_x = true_x + σ_loc * randn()
                obs_y = true_y + σ_loc * randn()

                # Random photon count
                photons = mean_photons * (0.5 + rand())
                σ_photons = sqrt(photons)

                # Create emitter with track_id = mol_id (ground truth)
                e = Emitter2DFit{Float64}(
                    obs_x, obs_y,
                    photons, 10.0,           # photons, bg
                    σ_loc, σ_loc, 0.0,       # σ_x, σ_y, σ_xy
                    σ_photons, 3.0,          # σ_photons, σ_bg
                    frame, 1, mol_id, 0      # frame, dataset, track_id, id
                )
                push!(emitters, e)
            end
        end
    end

    # Create camera and SMLD
    cam = IdealCamera(1:512, 1:512, 0.1)  # 0.1 μm pixels

    # Reset track_id to 0 for input (algorithm will populate it)
    input_emitters = [Emitter2DFit{Float64}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
        e.frame, e.dataset, 0, e.id  # track_id = 0
    ) for e in emitters]

    smld_input = BasicSMLD(input_emitters, cam, n_frames, 1)

    # Keep ground truth version
    smld_truth = BasicSMLD(emitters, cam, n_frames, 1)

    return smld_input, smld_truth, (mol_x, mol_y)
end

# Generate synthetic data
println("Generating synthetic SMLM data...")
smld_input, smld_truth, (true_x, true_y) = generate_synthetic_smld(
    n_molecules = 20,
    n_frames = 100,
    fov_size = 5.0,
    σ_loc = 0.02,
    p_on = 0.4
)

println("  Input localizations: $(length(smld_input.emitters))")
println("  True molecules: $(length(unique(e.track_id for e in smld_truth.emitters)))")

# Run frame connection
println("\nRunning frame connection...")
(combined, info) = frameconnect(
    smld_input;
    n_density_neighbors = 2,
    max_sigma_dist = 5.0,
    max_frame_gap = 5,
    max_neighbors = 2
)

println("  Combined localizations: $(length(combined.emitters))")
println("  Time: $(info.elapsed_s) seconds")

# Compare with ideal result (using ground truth track_id)
println("\nComputing ideal frame connection (ground truth)...")
smld_ideal_connected, smld_ideal_combined = defineidealFC(smld_truth; max_frame_gap = 5)
println("  Ideal combined: $(length(smld_ideal_combined.emitters))")

# Print estimated parameters
println("\nEstimated photophysics parameters:")
println("  k_on (dark→visible rate): $(round(info.k_on, digits=4)) /frame")
println("  k_off (visible→dark rate): $(round(info.k_off, digits=4)) /frame")
println("  k_bleach (bleaching rate): $(round(info.k_bleach, digits=4)) /frame")
println("  p_miss (miss probability): $(round(info.p_miss, digits=4))")

# Compute precision improvement
input_σ = mean([e.σ_x for e in smld_input.emitters])
combined_σ = mean([e.σ_x for e in combined.emitters])
println("\nPrecision improvement:")
println("  Mean input σ: $(round(input_σ * 1000, digits=1)) nm")
println("  Mean combined σ: $(round(combined_σ * 1000, digits=1)) nm")
println("  Improvement factor: $(round(input_σ / combined_σ, digits=2))x")

# Summary statistics
n_unique_tracks = length(unique(e.track_id for e in info.connected.emitters))
avg_locs_per_track = length(info.connected.emitters) / n_unique_tracks
println("\nConnection statistics:")
println("  Unique tracks identified: $n_unique_tracks")
println("  Avg localizations per track: $(round(avg_locs_per_track, digits=1))")
