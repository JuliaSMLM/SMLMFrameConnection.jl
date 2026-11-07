using Test, SMLMFrameConnection, SMLMData
include("test_helpers.jl")

using Random
using SMLMSim

# Same partition up to relabeling: i and j share a track in `a` iff they do in `b`.
function same_partition(a::Vector{Int}, b::Vector{Int})
    length(a) == length(b) || return false
    fwd = Dict{Int,Int}()
    bwd = Dict{Int,Int}()
    for (x, y) in zip(a, b)
        get!(fwd, x, y) == y || return false
        get!(bwd, y, x) == x || return false
    end
    return true
end

@testset "Input order" begin
    @testset "frameconnect does not depend on localization order" begin
        # SMLMSim returns localizations grouped by emitter, not by frame.
        Random.seed!(20261107)
        params = StaticSMLMConfig(density = 2.0, σ_psf = 0.13, minphotons = 50,
                                  ndatasets = 2, nframes = 400, framerate = 50.0)
        camera = IdealCamera(1:64, 1:64, 0.1)
        molecule = GenericFluor(photons = 3000.0, k_off = 3.0, k_on = 0.3)
        smld_native, _ = simulate(params; pattern = Nmer2D(n = 1, d = 0.0),
                                  molecule = molecule, camera = camera)
        # Distinct ids, so a combined localization's id names its source localization.
        em = [Emitter2DFit{Float64}(e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy,
                                    e.σ_photons, e.σ_bg, e.frame, e.dataset, e.track_id, i)
              for (i, e) in enumerate(smld_native.emitters)]
        n = length(em)
        @test n > 200
        @test !issorted([(e.dataset, e.frame) for e in em])   # the case under test

        withorder(p) = BasicSMLD(em[p], smld_native.camera, smld_native.n_frames,
                                 smld_native.n_datasets, copy(smld_native.metadata))
        p_sorted = sortperm(em; by = e -> (e.dataset, e.frame))
        orders = Dict(
            "native" => collect(1:n),
            "sorted" => p_sorted,
            "reversed" => reverse(p_sorted),
            "shuffled" => randperm(MersenneTwister(1), n),
        )
        cfg = FrameConnectConfig(calibration = CalibrationConfig())
        results = Dict(k => frameconnect(withorder(p), cfg) for (k, p) in orders)
        fields(c) = [(e.x, e.y, e.σ_x, e.σ_y, e.σ_xy, e.photons, e.frame, e.dataset, e.id,
                      e.track_id) for e in c.emitters]

        ref_combined, ref_info = results["sorted"]
        ref_tracks = zeros(Int, n)
        ref_tracks[orders["sorted"]] = [e.track_id for e in ref_info.connected.emitters]
        for (k, p) in orders
            combined, info = results[k]
            @testset "$k" begin
                @test info.n_combined == ref_info.n_combined
                @test info.n_preclusters == ref_info.n_preclusters
                @test info.k_on == ref_info.k_on
                @test info.k_off == ref_info.k_off
                @test info.k_bleach == ref_info.k_bleach
                @test info.p_miss == ref_info.p_miss
                # info.connected keeps the caller's order: map back to native order.
                tracks = zeros(Int, n)
                tracks[p] = [e.track_id for e in info.connected.emitters]
                @test same_partition(tracks, ref_tracks)
                @test [e.x for e in info.connected.emitters] == [e.x for e in em[p]]
                # The combined output, bit for bit, and the calibration.
                @test fields(combined) == fields(ref_combined)
                @test info.calibration.k_scale == ref_info.calibration.k_scale
                @test info.calibration.sigma_motion_nm == ref_info.calibration.sigma_motion_nm
            end
        end
    end

    @testset "combinelocalizations does not depend on input order" begin
        # Each combined localization takes frame and id from its track's earliest one.
        em = [make_emitter(5.0 + 0.001f, 5.0, f; track_id = 1) for f in (3, 1, 2)]
        em = [Emitter2DFit{Float64}(e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_xy,
                                    e.σ_photons, e.σ_bg, e.frame, e.dataset, e.track_id,
                                    10 * e.frame) for e in em]
        a = combinelocalizations(make_test_smld(em; n_frames = 3)).emitters[1]
        b = combinelocalizations(make_test_smld(reverse(em); n_frames = 3)).emitters[1]
        @test a.frame == 1 && a.id == 10
        @test (a.x, a.y, a.σ_x, a.frame, a.id) == (b.x, b.y, b.σ_x, b.frame, b.id)
    end

    @testset "precluster joins every cluster a localization links to" begin
        # sigma 10 nm, threshold 5 * (10 + 10) nm = 0.1 um. A (x 1.00, 0.98) and B
        # (x 1.15, 1.17) in frames 1-2; C at x 1.075 in frame 3 lies within 0.1 um of
        # A2 and B2. The preclusters are the connected components: one cluster.
        σ = 0.01
        emitters = [
            make_emitter(1.00, 5.0, 1; σ_pos = σ),   # A1
            make_emitter(1.15, 5.0, 1; σ_pos = σ),   # B1
            make_emitter(0.98, 5.0, 2; σ_pos = σ),   # A2
            make_emitter(1.17, 5.0, 2; σ_pos = σ),   # B2
            make_emitter(1.075, 5.0, 3; σ_pos = σ),  # C
        ]
        smld = make_test_smld(emitters; n_frames = 3)
        pre = SMLMFrameConnection.precluster(smld, SMLMFrameConnection.ParamStruct())
        ids = [e.track_id for e in pre.emitters]
        @test ids[2] == ids[4]          # B1 and B2 (0.02 um apart, adjacent frames)
        @test length(unique(ids)) == 1
    end
end
