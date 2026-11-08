using Test, SMLMFrameConnection, SMLMData
include("test_helpers.jl")

using Random
using SMLMSim
using Statistics

# Helper: deterministic pseudo-random from hash (no Random import needed)
_hashval(seed::Int, i::Int) = mod(hash(seed + i * 7919), UInt64) / typemax(UInt64)
# Box-Muller: two uniform -> one normal
function _hashnorm(seed::Int, i::Int)
    u1 = clamp(_hashval(seed, i), 1e-10, 1.0 - 1e-10)
    u2 = _hashval(seed, i + 1000)
    return sqrt(-2 * log(u1)) * cos(2π * u2)
end

@testset "Calibration" begin
    # Tracks of `nframes` consecutive localizations at fixed positions; photon counts
    # vary across tracks so the CRLB spans a range. Noise comes from SMLMSim's
    # apply_noise: the CRLB part with σ_psf (σ = σ_psf/√N), then, if σ_motion > 0, a
    # motion jitter of σ_motion per localization (drawn with photons = 1 so the draw
    # has exactly that σ). The σ_x, σ_y written to the data are the true CRLB divided
    # by k, so the data's CRLB is k times too small.
    function calibration_fixture(; k, σ_motion, ntracks = 2000, nframes = 6, seed = 1,
                                 photons = (20.0, 400.0))
        Random.seed!(seed)
        cam = make_test_camera()
        model = Emitter2DFit{Float64}[]
        for t in 1:ntracks
            x, y = 2.0 + 46.0 * rand(), 2.0 + 46.0 * rand()
            # Default: CRLB 29 nm down to 6.5 nm at σ_psf 0.13.
            N = photons[1] + (photons[2] - photons[1]) * rand()
            for f in 1:nframes
                push!(model, Emitter2DFit{Float64}(x, y, N, 10.0, 0.0, 0.0, 0.0, 1.0, 1.0,
                                                   f, 1, t, 0))
            end
        end
        smld = SMLMSim.apply_noise(BasicSMLD(model, cam, nframes, 1, Dict{String,Any}()), 0.13)
        em = smld.emitters
        if σ_motion > 0
            unit = [Emitter2DFit{Float64}(e.x, e.y, 1.0, e.bg, e.σ_x, e.σ_y, 0.0, 1.0, 1.0,
                                          e.frame, e.dataset, e.track_id, e.id) for e in em]
            jittered = SMLMSim.apply_noise(BasicSMLD(unit, cam, nframes, 1, Dict{String,Any}()),
                                           σ_motion).emitters
        else
            jittered = em
        end
        out = [Emitter2DFit{Float64}(j.x, j.y, e.photons, e.bg, e.σ_x / k, e.σ_y / k, 0.0,
                                     1.0, 1.0, e.frame, e.dataset, e.track_id, e.id)
               for (e, j) in zip(em, jittered)]
        return BasicSMLD(out, cam, nframes, 1, Dict{String,Any}())
    end

    @testset "CalibrationConfig construction" begin
        cc = CalibrationConfig()
        @test cc.clamp_k_to_one == true
        @test cc.filter_high_chi2 == false
        @test cc.chi2_filter_threshold == 6.0

        cc2 = CalibrationConfig(clamp_k_to_one=false, filter_high_chi2=true, chi2_filter_threshold=4.0)
        @test cc2.clamp_k_to_one == false
        @test cc2.filter_high_chi2 == true
        @test cc2.chi2_filter_threshold == 4.0
    end

    @testset "CalibrationConfig in FrameConnectConfig" begin
        config = FrameConnectConfig()
        @test config.calibration === nothing

        config2 = FrameConnectConfig(calibration=CalibrationConfig())
        @test config2.calibration isa CalibrationConfig
        @test config2.n_density_neighbors == 2
    end

    @testset "analyze_calibration with synthetic tracks" begin
        true_sigma_motion = 0.005  # 5 nm in μm
        true_k = 1.5
        smld = calibration_fixture(k = true_k, σ_motion = true_sigma_motion,
                                   ntracks = 200, nframes = 5)
        config = CalibrationConfig(clamp_k_to_one=false)
        result = SMLMFrameConnection.analyze_calibration(smld, config)

        @test result isa CalibrationResult
        @test result.calibration_applied == true
        @test result.warning == ""
        @test result.n_pairs > 0
        @test result.n_tracks_used > 0
        @test result.r_squared > 0.3

        @test abs(result.k_scale - true_k) < 0.25   # 800 pairs: sd of k is about 0.1
        @test result.sigma_motion_nm >= 0.0
        @test result.sigma_motion_nm < 20.0

        @test length(result.bin_centers) > 0
        @test length(result.bin_observed) == length(result.bin_centers)
        @test length(result.bin_counts) == length(result.bin_centers)
        @test all(result.bin_counts .>= 5)
        @test haskey(result.frame_shifts, 1)
    end

    @testset "analyze_calibration recovers k and σ_motion" begin
        # Before 2026-11 a fixed cut at chi² 6 on the uncalibrated CRLB gave k = 1.14
        # for a true 2.0, and the pair intercept 2σ_motion² was read as σ_motion²
        # (√2 too large once the cut is gone: 13.8 nm for a true 10).
        cfg = CalibrationConfig(clamp_k_to_one = false)
        for (k, σm_nm) in ((1.0, 0.0), (1.0, 10.0), (1.5, 5.0), (2.0, 0.0), (2.0, 5.0))
            r = SMLMFrameConnection.analyze_calibration(
                calibration_fixture(k = k, σ_motion = σm_nm / 1000), cfg)
            @testset "k = $k, σ_motion = $σm_nm nm" begin
                @test r.calibration_applied
                @test abs(r.k_scale - k) < 0.05 * k
                @test abs(r.sigma_motion_nm - σm_nm) < 2.5
                # The fit line goes through the stored (corrected) bins, within each
                # bin's noise: a bin mean of n pair values (dx²+dy²)/2 has relative
                # standard error about 1/√n.
                @test all(abs.(r.A .+ r.B .* r.bin_centers .- r.bin_observed) .<
                          5 .* r.bin_observed ./ sqrt.(r.bin_counts))
            end
        end
    end

    @testset "analyze_calibration does not apply an undetermined slope" begin
        # Precisions in a narrow range (CRLB 17.5-18.4 nm): intercept and slope are not
        # separable, and an unchecked fit landed on B < 0 with a 26 nm motion term.
        smld = calibration_fixture(k = 1.0, σ_motion = 0.0, photons = (50.0, 55.0))
        for clamp in (true, false)
            r = SMLMFrameConnection.analyze_calibration(smld,
                    CalibrationConfig(clamp_k_to_one = clamp))
            @test !r.calibration_applied ||
                  (abs(r.k_scale - 1) < 0.05 && r.sigma_motion_nm < 2.5)
            @test r.k_scale > 0
        end
    end

    @testset "applied variance matches the truth when a bound binds" begin
        # Per-localization variance the calibration applies at the median reported σ,
        # against the truth σ_motion² + k²·σ². A bound that binds (σ_motion² ≥ 0, or
        # k ≥ 1 with the clamp) refits the other parameter; keeping the jointly fitted
        # value instead gave 1.0-2.0 times the truth on these fixtures.
        cases = ((k = 1.0, σm = 0.0, photons = (50.0, 55.0)),
                 (k = 1.5, σm = 0.005, photons = (20.0, 400.0)))
        for c in cases, seed in 1:10, clamp in (true, false)
            smld = calibration_fixture(k = c.k, σ_motion = c.σm, photons = c.photons, seed = seed)
            r = SMLMFrameConnection.analyze_calibration(smld,
                    CalibrationConfig(clamp_k_to_one = clamp))
            σ = median(e.σ_x for e in smld.emitters)
            applied = r.calibration_applied ? (r.sigma_motion_nm / 1000)^2 + r.k_scale^2 * σ^2 : σ^2
            truth = c.σm^2 + c.k^2 * σ^2
            @test abs(applied / truth - 1) < 0.10
        end
    end

    @testset "analyze_calibration ignores mismatched pairs" begin
        # 1% of the pairs are wrong links far beyond the localization error.
        smld = calibration_fixture(k = 1.5, σ_motion = 0.0)
        em = copy(smld.emitters)
        Random.seed!(7)
        for i in 6:6:length(em)
            rand() < 0.06 || continue
            e = em[i]
            em[i] = Emitter2DFit{Float64}(e.x + 0.3, e.y, e.photons, e.bg, e.σ_x, e.σ_y,
                                          0.0, 1.0, 1.0, e.frame, e.dataset, e.track_id, e.id)
        end
        smld_bad = BasicSMLD(em, smld.camera, smld.n_frames, 1, Dict{String,Any}())
        r = SMLMFrameConnection.analyze_calibration(smld_bad,
                CalibrationConfig(clamp_k_to_one = false))
        @test r.calibration_applied
        @test abs(r.k_scale - 1.5) < 0.075
    end

    @testset "apply_calibration adds σ_motion² per localization" begin
        smld = calibration_fixture(k = 1.0, σ_motion = 0.010)
        r = SMLMFrameConnection.analyze_calibration(smld, CalibrationConfig(clamp_k_to_one = false))
        cal = SMLMFrameConnection.apply_calibration(smld, r)
        e0, e1 = smld.emitters[1], cal.emitters[1]
        @test e1.σ_x^2 ≈ (r.sigma_motion_nm / 1000)^2 + r.k_scale^2 * e0.σ_x^2
        @test r.A ≈ 2 * (r.sigma_motion_nm / 1000)^2 rtol = 1e-8
    end

    @testset "analyze_calibration with uniform brightness" begin
        # All emitters have same σ -- narrow CRLB range
        # With no motion and k=1, expect k_scale≈1
        emitters = Emitter2DFit{Float64}[]
        σ = 0.015
        n_tracks = 200
        seed = 9000

        for t in 1:n_tracks
            seed += 1
            x = 5.0 + _hashval(seed, 1) * 40.0
            y = 5.0 + _hashval(seed, 2) * 40.0
            for f in 1:4
                seed += 1
                # Jitter consistent with σ (no extra motion)
                dx = _hashnorm(seed, 1) * σ
                dy = _hashnorm(seed, 2) * σ
                push!(emitters, Emitter2DFit{Float64}(
                    x + dx, y + dy,
                    1000.0, 10.0, σ, σ, 0.0, 31.6, 3.16,
                    f, 1, t, 0
                ))
            end
        end

        smld = make_test_smld(emitters; n_frames=4)
        config = CalibrationConfig(clamp_k_to_one=false)
        result = SMLMFrameConnection.analyze_calibration(smld, config)

        # Should either succeed with reasonable k or fall back gracefully
        if result.calibration_applied
            @test result.k_scale > 0.5
            @test result.k_scale < 2.0
        else
            @test result.k_scale == 1.0  # fallback default
        end
    end

    @testset "analyze_calibration fallback on insufficient data" begin
        emitters = [
            make_emitter(5.0, 5.0, 1; track_id=1),
            make_emitter(5.01, 5.01, 2; track_id=1),
        ]
        smld = make_test_smld(emitters; n_frames=2)
        config = CalibrationConfig()
        result = SMLMFrameConnection.analyze_calibration(smld, config)

        @test result.calibration_applied == false
        @test !isempty(result.warning)
        @test result.k_scale == 1.0
    end

    @testset "apply_calibration modifies uncertainties" begin
        emitters = [
            make_emitter(5.0, 5.0, 1; σ_pos=0.02, track_id=1),
            make_emitter(5.01, 5.01, 2; σ_pos=0.03, track_id=1),
        ]
        smld = make_test_smld(emitters; n_frames=2)

        result = CalibrationResult(
            10.0, 1.5,
            2e-4, 2.25, 0.0, 0.0,
            0.95, 2.0,
            2, 1, 0,
            Float64[], Float64[], Int[],
            Dict{Int, Vector{NTuple{2,Float64}}}(),
            true, ""
        )

        smld_cal = SMLMFrameConnection.apply_calibration(smld, result)

        @test smld_cal.emitters[1].σ_x > smld.emitters[1].σ_x
        @test smld_cal.emitters[1].σ_y > smld.emitters[1].σ_y
        @test smld_cal.emitters[2].σ_x > smld.emitters[2].σ_x

        σ_motion = 10.0 / 1000
        k = 1.5
        expected_σx_1 = sqrt(σ_motion^2 + k^2 * 0.02^2)
        @test smld_cal.emitters[1].σ_x ≈ expected_σx_1 atol=1e-10

        @test smld_cal.emitters[1].x == smld.emitters[1].x
        @test smld_cal.emitters[1].y == smld.emitters[1].y
    end

    @testset "apply_calibration with sigma_xy" begin
        e = Emitter2DFit{Float64}(
            5.0, 5.0, 1000.0, 10.0,
            0.02, 0.03, 0.001,
            31.6, 3.16,
            1, 1, 1, 0
        )
        smld = make_test_smld([e]; n_frames=1)

        result = CalibrationResult(
            5.0, 2.0,
            5e-5, 4.0, 0.0, 0.0,
            0.95, 2.0,
            100, 50, 0,
            Float64[], Float64[], Int[],
            Dict{Int, Vector{NTuple{2,Float64}}}(),
            true, ""
        )

        smld_cal = SMLMFrameConnection.apply_calibration(smld, result)

        expected_σxy = 2.0^2 * 0.001
        @test smld_cal.emitters[1].σ_xy ≈ expected_σxy atol=1e-10
    end

    @testset "apply_calibration skips when not applied" begin
        emitters = [make_emitter(5.0, 5.0, 1; track_id=1)]
        smld = make_test_smld(emitters; n_frames=1)

        result = CalibrationResult(
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, NaN,
            0, 0, 0,
            Float64[], Float64[], Int[],
            Dict{Int, Vector{NTuple{2,Float64}}}(),
            false, "fallback"
        )

        smld_cal = SMLMFrameConnection.apply_calibration(smld, result)
        @test smld_cal === smld
    end

    @testset "frameconnect without calibration (regression)" begin
        smld = make_two_molecules_smld()
        (combined, info) = frameconnect(smld)

        @test info.calibration === nothing
        @test info.n_combined <= info.n_input
        @test length(combined.emitters) == info.n_combined
    end

    @testset "frameconnect with calibration" begin
        emitters = Emitter2DFit{Float64}[]
        n_molecules = 100
        frames_per = 4
        seed = 5000

        for m in 1:n_molecules
            seed += 1
            x = 5.0 + _hashval(seed, 1) * 40.0
            y = 5.0 + _hashval(seed, 2) * 40.0
            σ = 0.01 + _hashval(seed, 3) * 0.02

            for f in 1:frames_per
                seed += 1
                push!(emitters, Emitter2DFit{Float64}(
                    x + _hashnorm(seed, 1) * σ, y + _hashnorm(seed, 2) * σ,
                    1000.0, 10.0,
                    σ, σ, 0.0, 31.6, 3.16,
                    f, 1, 0, 0
                ))
            end
        end

        smld = make_test_smld(emitters; n_frames=frames_per)

        (combined_no_cal, info_no_cal) = frameconnect(smld)
        @test info_no_cal.calibration === nothing

        config = FrameConnectConfig(calibration=CalibrationConfig(clamp_k_to_one=false))
        (combined_cal, info_cal) = frameconnect(smld, config)

        @test info_cal.calibration isa CalibrationResult
        @test length(combined_no_cal.emitters) > 0
        @test length(combined_cal.emitters) > 0
    end

    @testset "chi2 track filtering" begin
        # Pre-assigned track_ids so we control exactly which tracks are outliers
        emitters = Emitter2DFit{Float64}[]

        # 150 normal tracks: positions barely move between frames, varying σ
        for t in 1:150
            x, y = 5.0 + t * 0.3, 5.0 + t * 0.2
            σ = 0.01 + 0.03 * (t / 150)  # σ varies from 0.01 to 0.04
            for f in 1:4
                # Tiny deterministic offset per frame
                dx = 0.001 * ((f % 3) - 1)
                push!(emitters, Emitter2DFit{Float64}(
                    x + dx, y + dx,
                    1000.0, 10.0, σ, σ, 0.0, 31.6, 3.16,
                    f, 1, t, 0
                ))
            end
        end

        # 10 outlier tracks: huge jump in frame 3
        for t in 1:10
            tid = 150 + t
            x, y = 5.0 + tid * 0.3, 5.0 + tid * 0.2
            σ = 0.015
            for f in 1:4
                jump = f >= 3 ? 5.0 : 0.0  # 5 μm jump = chi2 >> 6
                push!(emitters, Emitter2DFit{Float64}(
                    x + jump, y,
                    1000.0, 10.0, σ, σ, 0.0, 31.6, 3.16,
                    f, 1, tid, 0
                ))
            end
        end

        smld = make_test_smld(emitters; n_frames=4)

        config_no_filter = CalibrationConfig(filter_high_chi2=false)
        result_no_filter = SMLMFrameConnection.analyze_calibration(smld, config_no_filter)

        config_filter = CalibrationConfig(filter_high_chi2=true, chi2_filter_threshold=6.0)
        result_filter = SMLMFrameConnection.analyze_calibration(smld, config_filter)

        @test result_filter.n_tracks_filtered > 0
        @test result_filter.n_pairs < result_no_filter.n_pairs
    end
end
