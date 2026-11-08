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
    # by k, so the data's CRLB is k times too small. With ρ ≠ 0 the data are read in a
    # sheared coordinate frame, x' = x, y' = ρx + √(1-ρ²)y, applied to SMLMSim's
    # localizations before the motion jitter: the CRLB covariance then has correlation
    # ρ (σ_x' = σ_y' = σ, σ_xy' = ρσ²) and the motion stays isotropic.
    function calibration_fixture(; k, σ_motion, ntracks = 2000, nframes = 6, seed = 1,
                                 photons = (20.0, 400.0), ρ = 0.0)
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
        c = sqrt(1 - ρ^2)
        em = [Emitter2DFit{Float64}(e.x, ρ * e.x + c * e.y, e.photons, e.bg, e.σ_x, e.σ_y,
                                    ρ * e.σ_x * e.σ_y, 1.0, 1.0, e.frame, e.dataset,
                                    e.track_id, e.id) for e in smld.emitters]
        if σ_motion > 0
            unit = [Emitter2DFit{Float64}(e.x, e.y, 1.0, e.bg, e.σ_x, e.σ_y, 0.0, 1.0, 1.0,
                                          e.frame, e.dataset, e.track_id, e.id) for e in em]
            jittered = SMLMSim.apply_noise(BasicSMLD(unit, cam, nframes, 1, Dict{String,Any}()),
                                           σ_motion).emitters
        else
            jittered = em
        end
        out = [Emitter2DFit{Float64}(j.x, j.y, e.photons, e.bg, e.σ_x / k, e.σ_y / k,
                                     e.σ_xy / k^2, 1.0, 1.0, e.frame, e.dataset, e.track_id, e.id)
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

    @testset "applied variance matches the truth (bounds, narrow ranges)" begin
        # Per-localization variance the calibration applies at the median reported σ,
        # against the truth σ_motion² + k²·σ². A bound that binds (σ_motion² ≥ 0, or
        # k ≥ 1 with the clamp) refits the other parameter; keeping the jointly fitted
        # value instead gave 1.0-2.0 times the truth on these fixtures.
        cases = ((k = 1.0, σm = 0.0, photons = (50.0, 55.0)),
                 (k = 1.5, σm = 0.005, photons = (50.0, 55.0)),   # split not separable
                 (k = 1.5, σm = 0.005, photons = (20.0, 400.0)))
        for c in cases, seed in 1:10, clamp in (true, false)
            smld = calibration_fixture(k = c.k, σ_motion = c.σm, photons = c.photons, seed = seed)
            r = SMLMFrameConnection.analyze_calibration(smld,
                    CalibrationConfig(clamp_k_to_one = clamp))
            σ = median(e.σ_x for e in smld.emitters)
            applied = r.calibration_applied ? (r.sigma_motion_nm / 1000)^2 + r.k_scale^2 * σ^2 : σ^2
            truth = c.σm^2 + c.k^2 * σ^2
            @test r.calibration_applied
            @test abs(applied / truth - 1) < 0.10
        end
    end

    @testset "motion-dominated data get a motion term, not a scale" begin
        # Wide CRLB range (6.5-29 nm) with a dominant motion term: the observed variance
        # barely changes with the CRLB, so the free fit's slope is lost in the noise. A
        # scale-only fit there applied about half the variance to bright localizations
        # and 4-8 times to dim ones. Check the applied variance at the 2nd and 98th
        # percentile σ (or a fallback); the 800-pair case gets a looser bound (its
        # extremes reach about 23% over 40 seeds).
        for (σm, ntracks, tol) in ((0.100, 2000, 0.15), (0.050, 200, 0.25)), seed in 1:10
            smld = calibration_fixture(k = 1.0, σ_motion = σm, ntracks = ntracks, seed = seed)
            r = SMLMFrameConnection.analyze_calibration(smld, CalibrationConfig())
            r.calibration_applied || continue
            σs = sort([e.σ_x for e in smld.emitters])
            for σ in (σs[ceil(Int, 0.02 * length(σs))], σs[floor(Int, 0.98 * length(σs))])
                applied = (r.sigma_motion_nm / 1000)^2 + r.k_scale^2 * σ^2
                @test abs(applied / (σm^2 + σ^2) - 1) < tol
            end
        end
    end

    @testset "correlated CRLB covariances (σ_xy ≠ 0)" begin
        # The cut ignored σ_xy: its statistic was not chi²(2), more pairs fell past it
        # than the truncation factor assumes, and k came out low (1.85 for a true 2 at
        # ρ = 0.99; the uncorrelated control gives 1.99).
        cfg = CalibrationConfig(clamp_k_to_one = false)
        for (σm, tol) in ((0.0, 0.06), (0.004, 0.08))
            r = SMLMFrameConnection.analyze_calibration(
                    calibration_fixture(k = 2.0, σ_motion = σm, ρ = 0.99), cfg)
            @test r.calibration_applied
            @test abs(r.k_scale - 2.0) < tol
            # A trial fit with σ_motion² < 0 must not make the cut covariance
            # indefinite (that cut every pair and lost the free fit).
            @test !occursin("no free fit", r.warning)
        end
    end

    @testset "one reported precision: the one-parameter fits still apply" begin
        # Every localization has the same CRLB: no free fit exists (one bin), and the
        # calibration fell back to the identity with "Too few valid bins". Either
        # one-parameter fit gives the right variance at that σ.
        for (k, σm) in ((2.0, 0.0), (1.0, 0.010))
            smld = calibration_fixture(k = k, σ_motion = σm, photons = (100.0, 100.0))
            σ = smld.emitters[1].σ_x
            @test all(e -> e.σ_x == σ, smld.emitters)   # the case under test
            r = SMLMFrameConnection.analyze_calibration(smld,
                    CalibrationConfig(clamp_k_to_one = false))
            @test r.calibration_applied
            @test occursin("no free fit", r.warning)
            applied = (r.sigma_motion_nm / 1000)^2 + r.k_scale^2 * σ^2
            @test abs(applied / (σm^2 + k^2 * σ^2) - 1) < 0.03
        end
    end

    @testset "an applied calibration keeps at least 100 pairs" begin
        # Exactly 100 pairs: the free fit kept all 100, then a refit (bound or
        # one-parameter) ran its own cut and applied the calibration from 99 or fewer.
        cfg = CalibrationConfig(clamp_k_to_one = false)
        napplied = 0
        for seed in 1:30
            r = SMLMFrameConnection.analyze_calibration(
                    calibration_fixture(k = 1.5, σ_motion = 0.0, ntracks = 100, nframes = 2,
                                        seed = seed), cfg)
            r.calibration_applied || continue
            napplied += 1
            @test r.n_pairs >= 100
        end
        @test napplied > 0
    end

    @testset "R² is weighted about the weighted mean" begin
        # Three bins (reported 10, 20, 30; observed 1000, 1000, 1100; 100 pairs each).
        # The weighted R² used the unweighted mean of the bins: 0.7957 for 0.6231.
        reported = repeat([10.0, 20.0, 30.0]; inner = 100)
        observed = repeat([1000.0, 1000.0, 1100.0]; inner = 100)
        bin_of = repeat([1, 2, 3]; inner = 100)
        fit = SMLMFrameConnection._binned_wls(reported, observed, bin_of, trues(300), 3, 0.0, 1.0)
        c, o = [10.0, 20.0, 30.0], [1000.0, 1000.0, 1100.0]
        w = 1 ./ c .^ 2                      # model weights n/v², v = A + B·c = c
        X = hcat(ones(3), c)
        β = (X' * (w .* X)) \ (X' * (w .* o))
        r2 = 1 - sum(w .* (o .- X * β) .^ 2) / sum(w .* (o .- sum(w .* o) / sum(w)) .^ 2)
        @test fit.r_squared ≈ r2
        @test fit.r_squared ≈ 0.6231 atol = 1e-4
    end

    @testset "A and B report the free fit; the applied values are separate" begin
        # A clamped fit (true k 0.8 < 1): B is the fitted slope below 1, as before, so
        # consumers can label the clamp; B_applied is 1.
        smld = calibration_fixture(k = 0.8, σ_motion = 0.0)
        r = SMLMFrameConnection.analyze_calibration(smld, CalibrationConfig())
        @test r.calibration_applied
        @test r.B < 1
        @test abs(sqrt(r.B) - 0.8) < 0.05
        @test r.B_applied == 1.0 && r.k_scale == 1.0
        @test r.A_applied >= 0
        # Old 17-field constructor still works.
        r17 = CalibrationResult(2.0, 1.2, 1e-5, 1.44, 0.0, 0.0, 0.9, 2.0, 100, 10, 0,
                                Float64[], Float64[], Int[],
                                Dict{Int,Vector{NTuple{2,Float64}}}(), true, "")
        @test r17.B_applied ≈ 1.44 && r17.A_applied ≈ 2 * (0.002)^2
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
        @test r.A_applied ≈ 2 * (r.sigma_motion_nm / 1000)^2 rtol = 1e-8
        @test r.B_applied ≈ r.k_scale^2 rtol = 1e-12
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
            0.02, 0.03, 0.0003,           # covariance, < σ_x σ_y = 0.0006
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

        expected_σxy = 2.0^2 * 0.0003
        @test smld_cal.emitters[1].σ_xy ≈ expected_σxy atol=1e-10
        # apply_calibration returns on any input (a stepping loop never finished on
        # this non-positive-definite covariance, σ_xy 0.001 > σ_x σ_y 0.0006).
        bad = Emitter2DFit{Float64}(5.0, 5.0, 1000.0, 10.0, 0.02, 0.03, 0.001,
                                    31.6, 3.16, 1, 1, 1, 0)
        out = SMLMFrameConnection.apply_calibration(make_test_smld([bad]; n_frames=1), result)
        @test out.emitters[1].σ_xy ≈ 2.0^2 * 0.001 atol=1e-10
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

    @testset "the deprecated chi² track filter is ignored" begin
        # The filter dropped tracks with a pair chi² above the threshold under the
        # uncalibrated CRLB: the tracks that show k > 1. True k 1.5 gave k 1.05 at threshold
        # 6 (0.89 at 4) on this fixture.
        smld = calibration_fixture(k = 1.5, σ_motion = 0.0)
        off = SMLMFrameConnection.analyze_calibration(smld, CalibrationConfig())
        on = @test_logs (:warn, r"filter_high_chi2") match_mode = :any begin
            SMLMFrameConnection.analyze_calibration(smld,
                CalibrationConfig(filter_high_chi2 = true, chi2_filter_threshold = 6.0))
        end
        @test abs(on.k_scale - 1.5) < 0.05
        @test (on.k_scale, on.sigma_motion_nm, on.n_pairs) ==
              (off.k_scale, off.sigma_motion_nm, off.n_pairs)
        @test on.n_tracks_filtered == 0
    end
end
