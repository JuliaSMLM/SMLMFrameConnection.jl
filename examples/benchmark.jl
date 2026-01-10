# Performance benchmark for SMLMFrameConnection
#
# This script measures:
# - Execution time for different dataset sizes
# - Memory allocations
# - Scaling behavior
#
# Run with: julia --project examples/benchmark.jl

using SMLMFrameConnection
using SMLMData
using Random
using Statistics
using Printf

#=
Generate synthetic SMLM data for benchmarking.
Returns data with known ground truth for validation.
=#
function generate_benchmark_data(;
    n_molecules::Int,
    n_frames::Int,
    fov_size::Float64 = 10.0,
    σ_loc::Float64 = 0.02,
    p_on::Float64 = 0.3
)
    Random.seed!(12345)  # Fixed seed for reproducibility

    emitters = Emitter2DFit{Float64}[]

    mol_x = fov_size * rand(n_molecules)
    mol_y = fov_size * rand(n_molecules)

    for mol_id in 1:n_molecules
        true_x, true_y = mol_x[mol_id], mol_y[mol_id]

        for frame in 1:n_frames
            if rand() < p_on
                obs_x = true_x + σ_loc * randn()
                obs_y = true_y + σ_loc * randn()
                photons = 1000.0 * (0.5 + rand())

                e = Emitter2DFit{Float64}(
                    obs_x, obs_y,
                    photons, 10.0,
                    σ_loc, σ_loc,
                    sqrt(photons), 3.0,
                    frame, 1, 0, 0
                )
                push!(emitters, e)
            end
        end
    end

    cam = IdealCamera(1:512, 1:512, 0.1)
    return BasicSMLD(emitters, cam, n_frames, 1)
end

#=
Benchmark a single run with timing and allocation tracking.
=#
function benchmark_single(smld; warmup::Bool = false)
    # Warmup run (if needed)
    if warmup
        frameconnect(smld; nnearestclusters=1)
    end

    # Timed run
    GC.gc()  # Clean up before measurement

    stats = @timed frameconnect(smld; nnearestclusters=2, nsigmadev=5.0, maxframegap=5, nmaxnn=2)

    return (
        time = stats.time,
        bytes = stats.bytes,
        gctime = stats.gctime,
        n_input = length(smld.emitters),
        n_output = length(stats.value[3].emitters)
    )
end

#=
Run scaling benchmark across different dataset sizes.
=#
function run_scaling_benchmark()
    println("=" ^ 70)
    println("SMLMFrameConnection Performance Benchmark")
    println("=" ^ 70)

    # Test configurations: (n_molecules, n_frames)
    configs = [
        (10, 50),
        (20, 100),
        (50, 100),
        (100, 100),
        (100, 200),
        (200, 200),
    ]

    println("\nScaling Benchmark:")
    println("-" ^ 70)
    @printf("%-12s %-12s %-12s %-12s %-12s %-12s\n",
            "Molecules", "Frames", "Locs", "Time (s)", "Alloc (MB)", "Output")
    println("-" ^ 70)

    results = []

    for (n_mol, n_frames) in configs
        smld = generate_benchmark_data(n_molecules=n_mol, n_frames=n_frames)

        # Skip if too few localizations
        if length(smld.emitters) < 10
            println("Skipping ($n_mol, $n_frames) - too few localizations")
            continue
        end

        # Warmup on first run
        warmup = isempty(results)

        try
            result = benchmark_single(smld; warmup=warmup)

            @printf("%-12d %-12d %-12d %-12.3f %-12.1f %-12d\n",
                    n_mol, n_frames, result.n_input, result.time,
                    result.bytes / 1e6, result.n_output)

            push!(results, (
                n_molecules = n_mol,
                n_frames = n_frames,
                n_localizations = result.n_input,
                time = result.time,
                bytes = result.bytes,
                n_output = result.n_output
            ))
        catch e
            @printf("%-12d %-12d %-12s ERROR: %s\n",
                    n_mol, n_frames, "-", string(e)[1:min(30, length(string(e)))])
        end
    end

    println("-" ^ 70)

    return results
end

#=
Detailed profiling of a single run.
=#
function profile_detailed(; n_molecules::Int = 50, n_frames::Int = 100)
    println("\n" * "=" ^ 70)
    println("Detailed Allocation Profile")
    println("=" ^ 70)

    smld = generate_benchmark_data(n_molecules=n_molecules, n_frames=n_frames)
    println("\nDataset: $n_molecules molecules, $n_frames frames")
    println("Input localizations: $(length(smld.emitters))")

    # Warmup
    frameconnect(smld; nnearestclusters=1)

    println("\nRunning with @time macro:")
    println("-" ^ 70)

    GC.gc()
    @time smld_connected, smld_preclustered, smld_combined, params = frameconnect(
        smld;
        nnearestclusters = 2,
        nsigmadev = 5.0,
        maxframegap = 5,
        nmaxnn = 2
    )

    println("\nOutput:")
    println("  Preclusters: $(length(unique(e.track_id for e in smld_preclustered.emitters)))")
    println("  Final tracks: $(length(unique(e.track_id for e in smld_connected.emitters)))")
    println("  Combined emitters: $(length(smld_combined.emitters))")

    # Precision improvement
    input_σ = mean([e.σ_x for e in smld.emitters])
    combined_σ = mean([e.σ_x for e in smld_combined.emitters])
    println("\nPrecision improvement: $(round(input_σ / combined_σ, digits=2))x")

    return smld, smld_combined, params
end

#=
Memory allocation breakdown by function (approximate).
=#
function allocation_breakdown(; n_molecules::Int = 30, n_frames::Int = 50)
    println("\n" * "=" ^ 70)
    println("Allocation Breakdown (Approximate)")
    println("=" ^ 70)

    smld = generate_benchmark_data(n_molecules=n_molecules, n_frames=n_frames)
    println("\nDataset: $(length(smld.emitters)) localizations")

    # Warmup
    frameconnect(smld; nnearestclusters=1)

    # Prepare params
    params = SMLMFrameConnection.ParamStruct()
    params.nnearestclusters = 2
    params.nsigmadev = 5.0
    params.maxframegap = 5
    params.nmaxnn = 2

    println("\n@allocated for each step:")
    println("-" ^ 50)

    # Step 1: Precluster
    GC.gc()
    alloc_precluster = @allocated begin
        smld_pre = SMLMFrameConnection.precluster(smld, params)
    end
    smld_pre = SMLMFrameConnection.precluster(smld, params)
    @printf("  precluster:           %10.2f KB\n", alloc_precluster / 1024)

    # Step 2: Organize clusters
    GC.gc()
    alloc_organize = @allocated begin
        clusterdata = SMLMFrameConnection.organizeclusters(smld_pre)
    end
    clusterdata = SMLMFrameConnection.organizeclusters(smld_pre)
    @printf("  organizeclusters:     %10.2f KB\n", alloc_organize / 1024)

    # Step 3: Estimate params
    GC.gc()
    alloc_estimateparams = @allocated begin
        k_on, k_off, k_bleach, p_miss, _ = SMLMFrameConnection.estimateparams(smld_pre, clusterdata)
    end
    k_on, k_off, k_bleach, p_miss, _ = SMLMFrameConnection.estimateparams(smld_pre, clusterdata)
    params.k_on = k_on
    params.k_off = k_off
    params.k_bleach = k_bleach
    params.p_miss = p_miss
    @printf("  estimateparams:       %10.2f KB\n", alloc_estimateparams / 1024)

    # Step 4: Estimate densities
    GC.gc()
    alloc_densities = @allocated begin
        params.initialdensity = SMLMFrameConnection.estimatedensities(smld_pre, clusterdata, params)
    end
    params.initialdensity = SMLMFrameConnection.estimatedensities(smld_pre, clusterdata, params)
    @printf("  estimatedensities:    %10.2f KB\n", alloc_densities / 1024)

    # Step 5: Connect localizations (LAP)
    connectID_pre = [e.track_id for e in smld_pre.emitters]
    nframes = smld.n_frames > 0 ? smld.n_frames : maximum(e.frame for e in smld.emitters)
    GC.gc()
    alloc_connect = @allocated begin
        connectID = SMLMFrameConnection.connectlocalizations(connectID_pre, clusterdata, params, nframes)
    end
    @printf("  connectlocalizations: %10.2f KB\n", alloc_connect / 1024)

    # Step 6: Combine localizations
    # Need to create smld_connected first
    connectID = SMLMFrameConnection.connectlocalizations(connectID_pre, clusterdata, params, nframes)
    new_emitters = [SMLMData.Emitter2DFit{Float64}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg,
        e.frame, e.dataset, connectID[i], e.id
    ) for (i, e) in enumerate(smld.emitters)]
    smld_connected = BasicSMLD(new_emitters, smld.camera, smld.n_frames, smld.n_datasets)

    GC.gc()
    alloc_combine = @allocated begin
        smld_combined = combinelocalizations(smld_connected)
    end
    @printf("  combinelocalizations: %10.2f KB\n", alloc_combine / 1024)

    println("-" ^ 50)
    total = alloc_precluster + alloc_organize + alloc_estimateparams +
            alloc_densities + alloc_connect + alloc_combine
    @printf("  TOTAL:                %10.2f KB\n", total / 1024)
end

# Run all benchmarks
function main()
    println("\nSMLMFrameConnection.jl Performance Analysis")
    println("Julia version: $(VERSION)")
    println("Threads: $(Threads.nthreads())")
    println()

    # Scaling benchmark
    results = run_scaling_benchmark()

    # Detailed profile
    profile_detailed(n_molecules=50, n_frames=100)

    # Allocation breakdown
    allocation_breakdown(n_molecules=30, n_frames=50)

    println("\n" * "=" ^ 70)
    println("Benchmark complete!")
    println("=" ^ 70)

    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
