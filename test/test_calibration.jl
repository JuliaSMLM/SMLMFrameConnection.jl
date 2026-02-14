# Helper: deterministic pseudo-random from hash (no Random import needed)
_hashval(seed::Int, i::Int) = mod(hash(seed + i * 7919), UInt64) / typemax(UInt64)
# Box-Muller: two uniform -> one normal
function _hashnorm(seed::Int, i::Int)
    u1 = clamp(_hashval(seed, i), 1e-10, 1.0 - 1e-10)
    u2 = _hashval(seed, i + 1000)
    return sqrt(-2 * log(u1)) * cos(2π * u2)
end

@testset "Calibration" begin
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

        emitters = Emitter2DFit{Float64}[]
        track_id = 0
        n_tracks = 200
        frames_per_track = 5
        seed = 0

        for t in 1:n_tracks
            track_id += 1
            seed += 1
            base_x = 5.0 + _hashval(seed, 1) * 40.0
            base_y = 5.0 + _hashval(seed, 2) * 40.0
            crlb_sigma = 0.01 + _hashval(seed, 3) * 0.03

            for f in 1:frames_per_track
                seed += 1
                true_sigma = sqrt(true_sigma_motion^2 + true_k^2 * crlb_sigma^2)
                x = base_x + _hashnorm(seed, 1) * true_sigma
                y = base_y + _hashnorm(seed, 2) * true_sigma
                push!(emitters, Emitter2DFit{Float64}(
                    x, y, 1000.0, 10.0,
                    crlb_sigma, crlb_sigma, 0.0,
                    31.6, 3.16,
                    f, 1, track_id, 0
                ))
            end
        end

        smld = make_test_smld(emitters; n_frames=frames_per_track)
        config = CalibrationConfig(clamp_k_to_one=false)
        result = SMLMFrameConnection.analyze_calibration(smld, config)

        @test result isa CalibrationResult
        @test result.calibration_applied == true
        @test result.warning == ""
        @test result.n_pairs > 0
        @test result.n_tracks_used > 0
        @test result.r_squared > 0.3

        @test result.k_scale > 1.0
        @test result.k_scale < 3.0
        @test result.sigma_motion_nm >= 0.0
        @test result.sigma_motion_nm < 20.0

        @test length(result.bin_centers) > 0
        @test length(result.bin_observed) == length(result.bin_centers)
        @test haskey(result.frame_shifts, 1)
    end

    @testset "analyze_calibration with uniform brightness (ratio estimator)" begin
        # All emitters have same σ -- narrow CRLB range triggers ratio estimator
        # With no motion and k=1, expect k_scale≈1, sigma_motion≈0
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

        @test result.calibration_applied == true
        @test result.k_scale > 0.5
        @test result.k_scale < 2.0
        @test result.sigma_motion_nm == 0.0  # ratio estimator sets A=0
        @test result.A == 0.0
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
            Float64[], Float64[],
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
            Float64[], Float64[],
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
            Float64[], Float64[],
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
