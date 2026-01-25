# Realistic Frame Connection Example using SMLMSim
#
# This example demonstrates frame connection on a realistic-scale SMLM dataset
# generated with SMLMSim using physically-motivated blinking kinetics.
#
# Run with: julia --project=examples examples/realistic_frameconnection.jl

using SMLMFrameConnection
using SMLMSim
using SMLMData
using Statistics
using Printf

#=
Generate SMLM data with SMLMSim.

Key parameters from @sim:
- GenericFluor(; photons, k_off, k_on) for 2-state blinking
- k_off: on→off rate (Hz), controls blink duration
- k_on: off→on rate (Hz), controls dark time
- simulate() returns BasicSMLD with Emitter2DFit, track_id is ground truth
- Rates are per-second in SMLMSim, per-frame in frameconnect output

For realistic dSTORM:
- k_off = 2-5 Hz (blink duration)
- k_on = 0.01-0.05 Hz (dark time, gives 0.1-1% duty cycle)
- Density 10-100 emitters/μm²
=#

function main()
    println("=" ^ 70)
    println("Realistic Frame Connection Example (SMLMSim)")
    println("=" ^ 70)

    # =========================================================================
    # 1. Setup simulation parameters
    # =========================================================================
    println("\n1. Configuring SMLMSim parameters...")
    println("-" ^ 70)

    # Blinking kinetics (2-state model: on/off)
    # - k_off = 3 Hz → mean on-time ~333ms
    # - k_on = 0.03 Hz → mean off-time ~33s → ~1% duty cycle
    fluor = GenericFluor(;
        photons = 5000.0,    # Emission rate affects photon counts
        k_off = 3.0,         # On→off rate (Hz)
        k_on = 0.03          # Off→on rate (Hz) → 1% duty cycle
    )

    # Simulation parameters
    params = StaticSMLMParams(
        density = 20.0,       # Emitters per μm² (realistic dSTORM)
        σ_psf = 0.13,         # PSF width (μm) → precision = σ_psf/√photons
        minphotons = 100,     # Detection threshold
        nframes = 1000,       # Number of acquisition frames
        framerate = 50.0,     # Acquisition rate (Hz)
        ndatasets = 1
    )

    # Camera: 256×256 pixels at 100nm/pixel → 25.6×25.6 μm FOV
    camera = IdealCamera(1:256, 1:256, 0.1)

    @printf("  FOV size:        %.1f × %.1f μm\n", 25.6, 25.6)
    @printf("  Pixel size:      %.0f nm\n", 0.1 * 1000)
    @printf("  Density:         %.1f emitters/μm²\n", params.density)
    @printf("  Frames:          %d at %.0f fps\n", params.nframes, params.framerate)
    @printf("  PSF σ:           %.0f nm\n", params.σ_psf * 1000)

    println("\n  Blinking kinetics:")
    duty_cycle = 0.03 / (0.03 + 3.0) * 100
    @printf("    k_off:  %.1f Hz (mean on-time %.0f ms)\n", 3.0, 1000/3.0)
    @printf("    k_on:   %.2f Hz (mean off-time %.0f s)\n", 0.03, 1/0.03)
    @printf("    Duty cycle: %.1f%%\n", duty_cycle)

    # =========================================================================
    # 2. Run simulation
    # =========================================================================
    println("\n2. Running SMLMSim simulation...")
    println("-" ^ 70)

    smld_true, smld_model, smld_noisy = simulate(params;
        molecule = fluor,
        camera = camera,
        pattern = Nmer2D(n=1, d=0.0)  # Single isolated emitters
    )

    n_true_emitters = length(smld_true.emitters)
    n_localizations = length(smld_noisy.emitters)
    n_unique_tracks = length(unique(e.track_id for e in smld_noisy.emitters))

    @printf("  True emitters:      %d\n", n_true_emitters)
    @printf("  Localizations:      %d\n", n_localizations)
    @printf("  Unique tracks:      %d\n", n_unique_tracks)
    @printf("  Locs/emitter:       %.1f\n", n_localizations / n_true_emitters)

    # =========================================================================
    # 3. Prepare input for frame connection (clear track_id)
    # =========================================================================
    println("\n3. Preparing frame connection input...")
    println("-" ^ 70)

    # Clear track_id for input (algorithm should estimate it)
    input_emitters = [Emitter2DFit{Float64}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg,
        e.frame, e.dataset, 0, e.id  # track_id = 0
    ) for e in smld_noisy.emitters]

    smld_input = BasicSMLD(input_emitters, smld_noisy.camera,
                           smld_noisy.n_frames, smld_noisy.n_datasets)

    println("  track_id cleared for algorithm input")

    # =========================================================================
    # 4. Run frame connection
    # =========================================================================
    println("\n4. Running frame connection algorithm...")
    println("-" ^ 70)

    stats = @timed begin
        frameconnect(
            smld_input;
            nnearestclusters = 2,
            nsigmadev = 5.0,
            maxframegap = 10,  # Allow gaps since k_on=0.03 means reblinking
            nmaxnn = 2
        )
    end

    result = stats.value
    fc_params = result.params

    @printf("  Time:        %.2f seconds\n", stats.time)
    @printf("  Allocations: %.1f MB\n", stats.bytes / 1e6)

    # =========================================================================
    # 5. Results summary
    # =========================================================================
    n_tracks = length(unique(e.track_id for e in result.connected.emitters))
    n_combined = length(result.combined.emitters)

    println("\n5. Frame Connection Results")
    println("-" ^ 70)
    @printf("  Input localizations:  %d\n", length(smld_input.emitters))
    @printf("  Connected tracks:     %d\n", n_tracks)
    @printf("  Combined emitters:    %d\n", n_combined)
    @printf("  Reduction ratio:      %.1fx\n", length(smld_input.emitters) / n_combined)

    # =========================================================================
    # 6. Compare estimated vs true parameters
    # =========================================================================
    println("\n6. Parameter Estimation")
    println("-" ^ 70)

    # True per-frame rates (convert from Hz to per-frame at framerate)
    dt = 1.0 / params.framerate
    true_k_off_pf = 1 - exp(-3.0 * dt)    # k_off=3 Hz
    true_k_on_pf = 1 - exp(-0.03 * dt)    # k_on=0.03 Hz

    println("  Parameter    True (Hz)    Estimated (/frame)    True (/frame)")
    println("  " * "-" ^ 60)
    @printf("  k_on         %.2f         %.4f                %.6f\n",
            0.03, fc_params.k_on, true_k_on_pf)
    @printf("  k_off        %.2f         %.4f                %.4f\n",
            3.0, fc_params.k_off, true_k_off_pf)
    @printf("  k_bleach     N/A          %.5f                N/A (2-state model)\n",
            fc_params.k_bleach)
    @printf("  p_miss       -            %.4f                -\n", fc_params.p_miss)

    # =========================================================================
    # 7. Precision improvement
    # =========================================================================
    println("\n7. Precision Improvement")
    println("-" ^ 70)

    input_σ = mean([e.σ_x for e in smld_input.emitters])
    combined_σ = mean([e.σ_x for e in result.combined.emitters])
    avg_locs_per_track = length(result.connected.emitters) / n_tracks
    theoretical_improvement = sqrt(avg_locs_per_track)

    @printf("  Mean input σ:         %.1f nm\n", input_σ * 1000)
    @printf("  Mean combined σ:      %.1f nm\n", combined_σ * 1000)
    @printf("  Actual improvement:   %.2fx\n", input_σ / combined_σ)
    @printf("  Theoretical (√N):     %.2fx (N=%.1f locs/track)\n",
            theoretical_improvement, avg_locs_per_track)

    # =========================================================================
    # 8. Comparison with ground truth (ideal frame connection)
    # =========================================================================
    println("\n8. Ground Truth Comparison")
    println("-" ^ 70)

    # Use original smld_noisy which has true track_id from simulation
    smld_ideal_connected, smld_ideal_combined = defineidealFC(smld_noisy; maxframegap = 10)
    n_ideal = length(smld_ideal_combined.emitters)

    @printf("  True emitters (SMLMSim):    %d\n", n_true_emitters)
    @printf("  Ideal combined:             %d\n", n_ideal)
    @printf("  Algorithm combined:         %d\n", n_combined)

    if n_combined > n_ideal
        @printf("  Over-segmentation:          %.1f%% (algorithm found more)\n",
                100.0 * (n_combined - n_ideal) / n_ideal)
    else
        @printf("  Under-segmentation:         %.1f%% (algorithm merged too many)\n",
                100.0 * (n_ideal - n_combined) / n_ideal)
    end

    println("\n" * "=" ^ 70)
    println("Example complete!")
    println("=" ^ 70)

    return smld_input, result, smld_noisy
end

# Run
main()
