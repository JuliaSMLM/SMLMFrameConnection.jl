using Test, SMLMFrameConnection, SMLMData
include("test_helpers.jl")

@testset "Input checks and degenerate inputs" begin
    @testset "3D localizations are refused, not flattened" begin
        # x and y both vary, so the fixture is not also a zero-area precluster.
        em3 = [Emitter3DFit{Float64}(5.0 + 0.003i, 5.0 + 0.002i^2, 0.5, 1000.0, 10.0,
                                     0.01, 0.01, 0.03, 30.0, 3.0; frame = i)
               for i in 1:4]
        smld3 = BasicSMLD(em3, make_test_camera(), 4, 1, Dict{String,Any}())
        @test_throws ArgumentError frameconnect(smld3)
        @test_throws ArgumentError frameconnect(smld3, FrameConnectConfig())
        @test_throws ArgumentError combinelocalizations(smld3)
        @test_throws ArgumentError defineidealFC(smld3)
    end

    @testset "invalid inputs are refused at entry" begin
        good = [make_emitter(5.0 + 0.01i, 5.0 + 0.003i, i) for i in 1:4]
        bad(f) = (em = copy(good); em[2] = f(em[2]); make_test_smld(em; n_frames = 4))
        with(e; kw...) = Emitter2DFit{Float64}(get(kw, :x, e.x), e.y, e.photons, e.bg,
            get(kw, :σ_x, e.σ_x), e.σ_y, get(kw, :σ_xy, e.σ_xy), e.σ_photons, e.σ_bg,
            get(kw, :frame, e.frame), e.dataset, e.track_id, e.id)
        # Before: NaN σ crashed in precluster (BoundsError) or reached the LAP; zero σ
        # gave a SingularException or NaN combined positions; frame 0 and a too-small
        # n_frames gave a BoundsError in the cost matrix.
        for f in (e -> with(e; σ_x = NaN), e -> with(e; σ_x = 0.0), e -> with(e; x = NaN),
                  e -> with(e; σ_xy = 1.0), e -> with(e; frame = 0))
            @test_throws ArgumentError frameconnect(bad(f))
        end
        @test_throws ArgumentError frameconnect(make_test_smld(good; n_frames = 3))
        @test frameconnect(make_test_smld(good; n_frames = 0))[2].n_input == 4   # unknown
        # combinelocalizations: a zero σ gave NaN position and σ silently.
        tr = [make_emitter(5.0, 5.0, f; track_id = 1) for f in 1:2]
        tr[2] = with(tr[2]; σ_x = 0.0)
        @test_throws ArgumentError combinelocalizations(make_test_smld(tr; n_frames = 2))
    end

    @testset "combinelocalizations keeps unlinked and per-dataset tracks apart" begin
        # track_id 0 is SMLMData's "unlinked": three of them used to become one
        # localization at their mean position, across datasets.
        em = [make_emitter(1.0, 1.0, 1; dataset = 1), make_emitter(5.0, 5.0, 2; dataset = 1),
              make_emitter(9.0, 9.0, 3; dataset = 2)]
        out = combinelocalizations(make_test_smld(em; n_frames = 3, n_datasets = 2)).emitters
        @test length(out) == 3
        @test sort([e.x for e in out]) == [1.0, 5.0, 9.0]
        @test all(e -> e.track_id == 0, out)
        # One track_id in two datasets is two tracks.
        em = [make_emitter(1.0, 1.0, 1; track_id = 7, dataset = 1),
              make_emitter(1.0, 1.0, 2; track_id = 7, dataset = 1),
              make_emitter(9.0, 9.0, 1; track_id = 7, dataset = 2)]
        out = combinelocalizations(make_test_smld(em; n_frames = 2, n_datasets = 2)).emitters
        @test length(out) == 2
        @test sort([e.dataset for e in out]) == [1, 2]
    end

    @testset "one zero-area precluster returns" begin
        # Two localizations of one molecule at exactly one position: the precluster's
        # bounding box has zero area. This used to give an infinite density, a NaN
        # cost, and a LAP solve that never returned.
        smld = make_test_smld([make_emitter(5.0, 5.0, 1), make_emitter(5.0, 5.0, 2)];
                              n_frames = 2)
        combined, info = frameconnect(smld)
        @test info.n_combined == 1
        @test all(isfinite, info.initial_density)

        # Two same-frame localizations at one y (zero height): never connectable.
        smld = make_test_smld([make_emitter(5.0, 5.0, 1), make_emitter(5.05, 5.0, 1)];
                              n_frames = 1)
        combined, info = frameconnect(smld)
        @test info.n_combined == 2
    end

    # A cluster of two localizations two frames apart (frames 2 and 4, so their birth
    # costs are finite: the on-state density is 0 at frame 1), for create_costmatrix.
    function costcase(; initial_density, p_miss, frames = (2, 4))
        params = SMLMFrameConnection.ParamStruct()
        params.initial_density = [initial_density]
        params.k_on, params.k_off, params.k_bleach, params.p_miss = 0.01, 0.3, 0.01, p_miss
        # columns: x y σ_x σ_y σ_xy frame dataset connectID index
        cluster = [5.0 5.0 0.005 0.005 0.0 frames[1] 1 1 1;
                   5.001 5.0 0.005 0.005 0.0 frames[2] 1 1 2]
        k = params.k_on + params.k_off + params.k_bleach
        duty = params.k_on / k
        l1 = params.k_bleach * duty
        l2 = k - l1
        base = [duty * (exp(-l1 * (f - 1)) - exp(-l2 * (f - 1))) for f in 1:5]
        return SMLMFrameConnection.create_costmatrix([cluster], params, 1, 5, base)
    end

    @testset "non-finite costs never reach the LAP solver" begin
        # An infinite density times the zero on-state density at frame 1 gives NaN
        # in the birth block.
        cm = costcase(initial_density = Inf, p_miss = 0.1, frames = (1, 2))
        @test !any(isnan, cm)
    end

    @testset "a forbidden-by-p_miss link stays the most expensive" begin
        # p_miss = 0 makes the two-frame link cost Inf; it is replaced by a large value.
        # With a high density the birth costs are negative, so the replacement must not
        # be built from a signed sum (that made the forbidden link the cheapest).
        cm = costcase(initial_density = 1e12, p_miss = 0.0)
        link = cm[1, 2]
        finite_valid = [cm[2+1, 1], cm[2+2, 2], cm[1, 2+1], cm[2, 2+2]]   # birth, death
        @test isfinite(link)
        @test link > maximum(finite_valid)
        @test link > 0
        @test any(<(0), finite_valid)   # the case under test: negative costs present
    end

    @testset "integer columns of the cluster data stay exact above 2^24" begin
        # Track IDs (column 8) and input indices (column 9) are stored in the cluster
        # matrices; Float32 rounds integers above 2^24 = 16777216 to even neighbors,
        # which wrote track labels to the wrong localizations on large inputs.
        big = 2^24 + 1
        em = [make_emitter(5.0, 5.0, 1; track_id = big),
              make_emitter(9.0, 9.0, 1; track_id = big + 2)]
        smld = make_test_smld(em; n_frames = 1)
        clusters = SMLMFrameConnection.organizeclusters(smld)
        T = eltype(first(clusters))
        @test T(big) == big
        @test sort([Int(c[1, 8]) for c in clusters]) == [big, big + 2]
        @test sort([Int(c[1, 9]) for c in clusters]) == [1, 2]
    end
end
