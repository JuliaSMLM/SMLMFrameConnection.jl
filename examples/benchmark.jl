# Performance benchmark for SMLMFrameConnection
#
# Uses SMLMSim with SMLMAnalysis-equivalent parameters:
#   density=2.0/μm², 2000 frames, GenericFluor(photons=50000, k_off=20, k_on=0.04)
#   Camera: 256×128 px, 100 nm/px → 25.6×12.8 μm FOV
#
# Run with: julia --project=examples examples/benchmark.jl

using SMLMFrameConnection
using SMLMSim
using SMLMData
using Statistics
using Printf
using Random

"""
Generate a single dataset using SMLMSim with SMLMAnalysis-equivalent parameters.
Returns smld_noisy (no image generation).
"""
function generate_dataset(dataset_id::Int;
    nframes::Int = 2000,
    ndatasets::Int = 1,
    density::Float64 = 2.0,
    nx::Int = 256,
    ny::Int = 128,
    pixel_size::Float64 = 0.1,
    photons::Float64 = 50000.0,
    k_off::Float64 = 20.0,
    k_on::Float64 = 0.04,
    σ_psf::Float64 = 0.13
)
    Random.seed!(dataset_id * 42)

    fluor = GenericFluor(; photons=photons, k_off=k_off, k_on=k_on)
    camera = IdealCamera(nx, ny, pixel_size)

    sim_params = StaticSMLMConfig(;
        density=density,
        σ_psf=σ_psf,
        minphotons=100,
        nframes=nframes,
        ndatasets=ndatasets
    )

    smld_noisy, _ = simulate(sim_params;
        molecule=fluor,
        camera=camera,
        pattern=Nmer2D(n=1, d=0.0)
    )

    # Set dataset field and clear track_id
    new_emitters = [Emitter2DFit{Float64}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
        e.frame, dataset_id, 0, e.id
    ) for e in smld_noisy.emitters]

    return BasicSMLD(new_emitters, camera, nframes, 1)
end

function bench(f, n=10)
    f()  # warmup
    GC.gc()
    times = Float64[]
    allocs = Int[]
    for _ in 1:n
        s = @timed f()
        push!(times, s.time)
        push!(allocs, s.bytes)
    end
    return (median=median(times), min=minimum(times), mean=mean(times), alloc_bytes=median(allocs))
end

function cluster_stats(smld, params)
    smld_pre = SMLMFrameConnection.precluster(smld, params)
    clusterdata = SMLMFrameConnection.organizeclusters(smld_pre)
    sizes = [size(c, 1) for c in clusterdata]
    multi = filter(s -> s > 1, sizes)
    return (
        n_clusters = length(sizes),
        n_multi = length(multi),
        max_size = isempty(multi) ? 0 : maximum(multi),
        median_size = isempty(multi) ? 0.0 : median(multi),
        mean_size = isempty(multi) ? 0.0 : mean(multi)
    )
end

function profile_stages(smld)
    params = SMLMFrameConnection.ParamStruct()
    params.n_density_neighbors = 2; params.max_sigma_dist = 5.0
    params.max_frame_gap = 5; params.max_neighbors = 2

    # Cluster stats
    cs = cluster_stats(smld, params)
    @printf("  Clusters: %d total, %d multi-loc\n", cs.n_clusters, cs.n_multi)
    @printf("  Multi-loc: max=%d, median=%.1f, mean=%.1f\n", cs.max_size, cs.median_size, cs.mean_size)

    # Prepare intermediates for stage benchmarks
    smld_pre = SMLMFrameConnection.precluster(smld, params)
    clusterdata = SMLMFrameConnection.organizeclusters(smld_pre)
    k_on, k_off, k_bleach, p_miss, _ = SMLMFrameConnection.estimateparams(smld_pre, clusterdata)
    params.k_on = k_on; params.k_off = k_off; params.k_bleach = k_bleach; params.p_miss = p_miss
    params.initial_density = SMLMFrameConnection.estimatedensities(smld_pre, clusterdata, params)
    connectID_pre = [e.track_id for e in smld_pre.emitters]
    nframes = smld.n_frames
    connectID = SMLMFrameConnection.connectlocalizations(connectID_pre, clusterdata, params, nframes)
    new_emitters = [SMLMData.Emitter2DFit{Float64}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
        e.frame, e.dataset, connectID[i], e.id
    ) for (i, e) in enumerate(smld.emitters)]
    smld_connected = BasicSMLD(new_emitters, smld.camera, smld.n_frames, smld.n_datasets)

    stages = [
        ("precluster",           () -> SMLMFrameConnection.precluster(smld, params)),
        ("organizeclusters",     () -> SMLMFrameConnection.organizeclusters(smld_pre)),
        ("estimateparams",       () -> SMLMFrameConnection.estimateparams(smld_pre, clusterdata)),
        ("estimatedensities",    () -> SMLMFrameConnection.estimatedensities(smld_pre, clusterdata, params)),
        ("connectlocalizations", () -> SMLMFrameConnection.connectlocalizations(connectID_pre, clusterdata, params, nframes)),
        ("combinelocalizations", () -> combinelocalizations(smld_connected)),
    ]

    results = []
    for (name, f) in stages
        r = bench(f, 10)
        push!(results, (name=name, r...))
    end

    total_ms = sum(r.median for r in results) * 1000
    println()
    for r in results
        t_ms = r.median * 1000
        pct = 100 * t_ms / total_ms
        @printf("  %-25s  %8.3f ms  (%5.1f%%)  alloc: %8.1f KB\n",
                r.name, t_ms, pct, r.alloc_bytes / 1024)
    end
    println("  " * "-" ^ 60)
    @printf("  %-25s  %8.3f ms\n", "TOTAL (stages)", total_ms)

    # Full pipeline
    rf = bench(() -> frameconnect(smld; max_frame_gap=5, calibration=CalibrationConfig()), 10)
    @printf("  %-25s  %8.3f ms  alloc: %8.1f KB\n\n",
            "frameconnect() full", rf.median * 1000, rf.alloc_bytes / 1024)

    return results
end

function main()
    println("SMLMFrameConnection Benchmark (SMLMSim-generated data)")
    println("Julia $(VERSION), $(Threads.nthreads()) threads")
    println()

    # Warmup with small dataset
    print("Warmup... ")
    warmup = generate_dataset(0; nframes=200, density=0.5, nx=64, ny=64)
    frameconnect(warmup; max_frame_gap=5)
    println("done.\n")

    # Scaling configs: (nframes, ndatasets_to_combine, density)
    configs = [
        (2000, 1,  2.0,  "1×2k frames, ρ=2"),
        (2000, 4,  2.0,  "4×2k frames, ρ=2"),
        (5000, 1,  2.0,  "1×5k frames, ρ=2"),
        (2000, 1,  5.0,  "1×2k frames, ρ=5"),
        (5000, 4,  2.0,  "4×5k frames, ρ=2"),
    ]

    for (nf, nds, dens, label) in configs
        # Generate and optionally merge datasets
        smlds = [generate_dataset(ds; nframes=nf, density=dens) for ds in 1:nds]
        if nds == 1
            smld = smlds[1]
        else
            all_emitters = vcat([s.emitters for s in smlds]...)
            smld = BasicSMLD(all_emitters, smlds[1].camera, nf, nds)
        end

        println("=" ^ 70)
        @printf("%-30s → %d localizations\n", label, length(smld.emitters))
        println("=" ^ 70)
        profile_stages(smld)
    end

    println("Benchmark complete.")
end

main()
