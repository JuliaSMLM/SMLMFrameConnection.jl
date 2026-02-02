# Full-Scale Frame Connection Benchmark
#
# Stress-tests frame connection on realistic full-scale dSTORM data:
# - 20 datasets, 5000 frames each
# - 30 emitters/μm² density (~20,000 emitters per dataset)
# - ~0.5% duty cycle (realistic dSTORM blinking)
#
# Run with: julia --project=examples examples/fullscale_benchmark.jl

using SMLMFrameConnection
using SMLMSim
using SMLMData
using Statistics
using Printf
using Random

#=
Realistic dSTORM benchmark configuration:
- FOV: 25.6 × 25.6 μm (256×256 at 100nm/px)
- Density: 30/μm² (~20,000 emitters per dataset)
- Framerate: 50 Hz (typical dSTORM)
- Blinking: k_off=3 Hz, k_on=0.02 Hz → 0.66% duty cycle
- Expected: ~66k localizations per dataset (density × area × frames × duty_cycle)
=#

"""
    generate_dataset(dataset_id; kwargs...) -> BasicSMLD

Generate a single dataset using SMLMSim with realistic dSTORM parameters.
"""
function generate_dataset(dataset_id::Int;
    n_frames::Int = 5000,
    framerate::Float64 = 50.0,
    density::Float64 = 30.0,
    fov_pixels::Int = 256,
    pixel_size::Float64 = 0.1,
    k_off::Float64 = 3.0,
    k_on::Float64 = 0.02,
    photons::Float64 = 5000.0,
    σ_psf::Float64 = 0.13
)
    # Set seed for reproducibility per dataset
    Random.seed!(dataset_id * 12345)

    fluor = GenericFluor(; photons=photons, k_off=k_off, k_on=k_on)

    params = StaticSMLMParams(
        density = density,
        σ_psf = σ_psf,
        minphotons = 100,
        nframes = n_frames,
        framerate = framerate,
        ndatasets = 1
    )

    camera = IdealCamera(1:fov_pixels, 1:fov_pixels, pixel_size)

    _, _, smld_noisy = simulate(params;
        molecule = fluor,
        camera = camera,
        pattern = Nmer2D(n=1, d=0.0)
    )

    # Update dataset field to match dataset_id
    new_emitters = [Emitter2DFit{Float64}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg,
        e.frame, dataset_id, 0, e.id  # Clear track_id, set dataset
    ) for e in smld_noisy.emitters]

    return BasicSMLD(new_emitters, camera, n_frames, 1)
end

"""
    benchmark_single_dataset(smld) -> NamedTuple

Run frame connection on a single dataset and return timing/results.
"""
function benchmark_single_dataset(smld::BasicSMLD)
    n_input = length(smld.emitters)

    GC.gc()  # Clean before measurement

    stats = @timed begin
        frameconnect(
            smld;
            nnearestclusters = 2,
            nsigmadev = 5.0,
            maxframegap = 10,
            nmaxnn = 2
        )
    end

    (combined, info) = stats.value

    return (
        n_input = n_input,
        n_output = info.n_combined,
        n_tracks = info.n_tracks,
        time = stats.time,
        bytes = stats.bytes,
        gctime = stats.gctime
    )
end

"""
    run_fullscale_benchmark(; n_datasets=20, kwargs...)

Run the complete full-scale benchmark.
"""
function run_fullscale_benchmark(;
    n_datasets::Int = 20,
    n_frames::Int = 5000,
    density::Float64 = 30.0,
    warmup::Bool = true
)
    println("=" ^ 75)
    println("FULL-SCALE FRAME CONNECTION BENCHMARK")
    println("=" ^ 75)

    # Configuration
    fov_size = 25.6  # μm (256 pixels at 100nm)
    fov_area = fov_size^2
    expected_emitters = round(Int, density * fov_area)

    println("\nConfiguration:")
    println("-" ^ 75)
    @printf("  Datasets:          %d\n", n_datasets)
    @printf("  Frames/dataset:    %d\n", n_frames)
    @printf("  FOV:               %.1f × %.1f μm (256×256 px)\n", fov_size, fov_size)
    @printf("  Density:           %.1f emitters/μm²\n", density)
    @printf("  Expected emitters: ~%d per dataset\n", expected_emitters)
    @printf("  Blinking:          k_off=3 Hz, k_on=0.02 Hz (0.66%% duty cycle)\n")
    @printf("  Framerate:         50 Hz\n")

    # Warmup run (small dataset to compile)
    if warmup
        println("\nWarmup run (compiling)...")
        warmup_smld = generate_dataset(0; n_frames=500, density=1.0)
        println("  Generated $(length(warmup_smld.emitters)) warmup locs")
        if length(warmup_smld.emitters) > 10
            _ = benchmark_single_dataset(warmup_smld)
        end
        GC.gc()
        println("  Warmup complete.")
    end

    # Main benchmark
    println("\n" * "=" ^ 75)
    println("BENCHMARK RESULTS")
    println("=" ^ 75)
    println()
    @printf("%-8s  %10s  %10s  %10s  %10s  %12s\n",
            "Dataset", "Input", "Output", "Tracks", "Time (s)", "Alloc (MB)")
    println("-" ^ 75)

    results = []
    total_input = 0
    total_output = 0
    total_time = 0.0
    total_bytes = 0

    for ds in 1:n_datasets
        # Generate dataset
        print("  Generating dataset $ds... ")
        gen_stats = @timed generate_dataset(ds; n_frames=n_frames, density=density)
        smld = gen_stats.value
        @printf("(%d locs, %.1fs)\n", length(smld.emitters), gen_stats.time)

        # Benchmark frame connection
        result = benchmark_single_dataset(smld)
        push!(results, result)

        total_input += result.n_input
        total_output += result.n_output
        total_time += result.time
        total_bytes += result.bytes

        @printf("%-8d  %10d  %10d  %10d  %10.2f  %12.1f\n",
                ds, result.n_input, result.n_output, result.n_tracks,
                result.time, result.bytes / 1e6)

        # Free memory
        smld = nothing
        GC.gc()
    end

    # Summary statistics
    println("-" ^ 75)
    @printf("%-8s  %10d  %10d  %10s  %10.2f  %12.1f\n",
            "TOTAL", total_input, total_output, "-", total_time, total_bytes / 1e6)

    # Detailed summary
    times = [r.time for r in results]
    inputs = [r.n_input for r in results]
    outputs = [r.n_output for r in results]

    println("\n" * "=" ^ 75)
    println("SUMMARY STATISTICS")
    println("=" ^ 75)

    println("\nScale:")
    @printf("  Total localizations:     %d (%.2fM)\n", total_input, total_input/1e6)
    @printf("  Total combined:          %d\n", total_output)
    @printf("  Overall reduction:       %.1fx\n", total_input / total_output)
    @printf("  Avg locs per dataset:    %d\n", round(Int, mean(inputs)))

    println("\nTiming:")
    @printf("  Total time:              %.2f seconds\n", total_time)
    @printf("  Avg time per dataset:    %.2f seconds\n", mean(times))
    @printf("  Min/Max time:            %.2f / %.2f seconds\n", minimum(times), maximum(times))
    @printf("  Throughput:              %.0f locs/second\n", total_input / total_time)

    println("\nMemory:")
    @printf("  Total allocated:         %.1f GB\n", total_bytes / 1e9)
    @printf("  Avg per dataset:         %.1f MB\n", mean([r.bytes for r in results]) / 1e6)

    println("\nPer-localization metrics:")
    @printf("  Time per 1k locs:        %.2f ms\n", 1000 * total_time / (total_input / 1000))
    @printf("  Memory per 1k locs:      %.1f MB\n", total_bytes / 1e6 / (total_input / 1000))

    println("\n" * "=" ^ 75)
    println("Benchmark complete!")
    println("=" ^ 75)

    return results
end

# Run benchmark
results = run_fullscale_benchmark(
    n_datasets = 20,
    n_frames = 5000,
    density = 30.0,
    warmup = true
)
