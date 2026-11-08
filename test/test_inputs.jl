using Test, SMLMFrameConnection, SMLMData
include("test_helpers.jl")

using Random

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

    @testset "a positive-definite Float32 covariance is accepted" begin
        # σ_xy just below σ_xσ_y: positive definite (exact determinant 3.5e-15), but
        # σ_xy² and σ_x²σ_y² rounded equal in Float32 and the entry check refused it.
        sx, sy, sxy = 0.03f0, 0.02f0, 0.00059999997f0
        @test Float64(sxy) < Float64(sx) * Float64(sy)          # the case under test
        em = [Emitter2DFit{Float32}(5f0 + 0.002f0 * i, 5f0, 1000f0, 10f0, sx, sy, sxy,
                                    30f0, 3f0, i, 1, 0, 0) for i in 1:4]
        smld = BasicSMLD(em, make_test_camera(), 4, 1, Dict{String,Any}())
        combined, info = frameconnect(smld)
        @test info.n_input == 4
        @test all(e -> isfinite(e.x) && isfinite(e.σ_x), combined.emitters)
        # Singular in exact arithmetic (σ_xy = σ_xσ_y exactly) is still refused.
        sing = [Emitter2DFit{Float32}(5f0 + 0.002f0 * i, 5f0, 1000f0, 10f0, 0.5f0, 0.25f0,
                                      0.125f0, 30f0, 3f0, i, 1, 0, 0) for i in 1:2]
        @test_throws ArgumentError frameconnect(BasicSMLD(sing, make_test_camera(), 2, 1,
                                                          Dict{String,Any}()))
    end

    @testset "apply_calibration keeps covariances positive definite" begin
        # A correlation as close to 1 as Float32 holds, scaled by k² and rounded back to
        # Float32, crossed 1 in about 7% of cases; combinelocalizations then refused
        # the calibrated data.
        Random.seed!(11)
        em = Emitter2DFit{Float32}[]
        for i in 1:2000
            sx, sy = Float32(0.005 + 0.03rand()), Float32(0.005 + 0.03rand())
            sxy = prevfloat(Float32(Float64(sx) * Float64(sy)))
            while Float64(sxy) >= Float64(sx) * Float64(sy)
                sxy = prevfloat(sxy)
            end
            push!(em, Emitter2DFit{Float32}(5f0, 5f0, 1000f0, 10f0, sx, sy, sxy, 30f0, 3f0,
                                            1, 1, i, i))
        end
        smld = BasicSMLD(em, make_test_camera(), 1, 1, Dict{String,Any}())
        for ksq in (1.21, 1.4217910694289009, 2.0)
            r = CalibrationResult(0.0, sqrt(ksq), 0.0, ksq, 0.0, 0.0, 1.0, 2.0, 1000, 100, 0,
                                  Float64[], Float64[], Int[],
                                  Dict{Int,Vector{NTuple{2,Float64}}}(), true, "")
            cal = SMLMFrameConnection.apply_calibration(smld, r)
            @test all(e -> abs(Float64(e.σ_xy)) < Float64(e.σ_x) * Float64(e.σ_y), cal.emitters)
        end
    end

    @testset "a group of one keeps its covariance; small σ is not floored" begin
        # σ = 0.1 nm: σ_x²σ_y² = 1e-16 μm⁴ fell below an absolute floor eps() = 2.2e-16,
        # so a lone localization came back with σ 0.149 nm and a pair of them with
        # 0.105 nm instead of 0.0707 nm.
        em = [make_emitter(5.0, 5.0, 1; σ_pos = 1e-4),
              make_emitter(6.0, 6.0, 1; σ_pos = 1e-4, track_id = 1),
              make_emitter(6.0, 6.0, 2; σ_pos = 1e-4, track_id = 1)]
        out = sort(combinelocalizations(make_test_smld(em; n_frames = 2)).emitters;
                   by = e -> e.track_id)
        @test out[1].σ_x == 1e-4 && out[1].σ_y == 1e-4 && out[1].σ_xy == 0.0
        @test out[2].σ_x ≈ 1e-4 / sqrt(2) rtol = 1e-12
        @test out[2].σ_y ≈ 1e-4 / sqrt(2) rtol = 1e-12
        # A lone Float32 localization comes back bit for bit.
        e32 = Emitter2DFit{Float32}(5.1f0, 4.9f0, 1000f0, 10f0, 0.013f0, 0.011f0, 2f-5,
                                    30f0, 3f0, 1, 1, 0, 0)
        o32 = combinelocalizations(BasicSMLD([e32], make_test_camera(), 1, 1,
                                             Dict{String,Any}())).emitters[1]
        @test (o32.x, o32.y, o32.σ_x, o32.σ_y, o32.σ_xy) == (e32.x, e32.y, e32.σ_x, e32.σ_y, e32.σ_xy)
    end

    @testset "nearly singular covariances stay finite and positive definite" begin
        # Positive definite in exact arithmetic, |σ_xy| within a rounding error of
        # σ_xσ_y (inputs from the Codex re-review of cba16d3, where each failed).
        e64(x, y, f, t; tid = 0) = Emitter2DFit{Float64}(x, y, 1000.0, 10.0, t..., 30.0, 3.0,
                                                         f, 1, tid, 0)
        e32(x, y, f, t; tid = 0) = Emitter2DFit{Float32}(x, y, 1000f0, 10f0, t..., 30f0, 3f0,
                                                         f, 1, tid, 0)
        posdef(e) = SMLMFrameConnection._cov_det(e.σ_x, e.σ_y, e.σ_xy) > 0
        # Float64 σ_xσ_y rounds onto σ_xy: refused before, accepted now.
        t6 = (0.03, 0.02, 0.0006)
        @test SMLMFrameConnection._cov_det(t6...) > 0
        @test frameconnect(make_test_smld([e64(5.0 + 0.002i, 5.0, i, t6) for i in 1:4];
                                          n_frames = 4))[2].n_input == 4
        # The sum of two such covariances: create_costmatrix took log of a negative
        # determinant (DomainError).
        t4 = (0.025553879796731797, 0.018355861103728575, 0.00046906346821018463)
        p = SMLMFrameConnection.ParamStruct(); p.initial_density = [1.0]
        p.k_on, p.k_off, p.k_bleach, p.p_miss = 0.01, 0.3, 0.01, 0.1
        cl = [5.0 5.0 t4[1] t4[2] t4[3] 2 1 1 1; 5.001 5.0 t4[1] t4[2] t4[3] 3 1 1 2]
        cm = SMLMFrameConnection.create_costmatrix([cl], p, 1, 5, fill(0.03, 5))
        @test isfinite(cm[1, 2])
        # Combining two at one position gave (NaN, NaN, Inf, Inf, Inf).
        tb = (0.01761696674838718, 0.011928085556888963, 0.00021013668662763022)
        o = combinelocalizations(make_test_smld([e64(5.0, 5.0, 1, tb; tid = 1),
                                                 e64(5.0, 5.0, 2, tb; tid = 1)];
                                                n_frames = 2)).emitters[1]
        @test (o.x, o.y) == (5.0, 5.0)
        @test all(isfinite, (o.σ_x, o.σ_y, o.σ_xy)) && posdef(o)
        # Float32 rounding of the combined covariance crossed ρ = 1, so the output was
        # refused when combined again.
        t5 = (0.023896482f0, 0.028789306f0, 0.0006879631f0)
        o = combinelocalizations(make_test_smld([e32(5f0, 5f0, 1, t5; tid = 1),
                                                 e32(5f0, 5f0, 2, t5; tid = 1)];
                                                n_frames = 2)).emitters[1]
        @test posdef(o)
        @test length(combinelocalizations(make_test_smld([o]; n_frames = 2)).emitters) == 1
        # σ² underflowing in Float32 made apply_calibration loop forever.
        tiny = BasicSMLD([e32(5f0, 5f0, 1, (1f-24, 1f-24, 0f0))], make_test_camera(), 1, 1,
                         Dict{String,Any}())
        r = CalibrationResult(0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 1000, 100, 0,
                              Float64[], Float64[], Int[],
                              Dict{Int,Vector{NTuple{2,Float64}}}(), true, "")
        out = SMLMFrameConnection.apply_calibration(tiny, r).emitters[1]
        @test out.σ_x == 1f-24 && posdef(out)
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
