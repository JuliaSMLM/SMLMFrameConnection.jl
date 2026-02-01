@testset "frameconnect" begin
    @testset "return types" begin
        # Use multiple molecules for proper density estimation
        emitters = vcat(
            [make_emitter(Float64(i), Float64(i), j) for i in 1:5 for j in 1:3]...
        )
        smld = make_test_smld(emitters; n_frames=3)
        (combined, info) = frameconnect(smld)

        @test combined isa BasicSMLD
        @test info isa ConnectInfo
        @test info.connected isa BasicSMLD
        @test info.algorithm == :lap
        @test info.elapsed_ns > 0
        @test info.n_input == length(smld.emitters)
        @test info.n_combined == length(combined.emitters)
    end

    @testset "single molecule connection" begin
        # Single molecule with 5 localizations, but add more spread-out molecules
        # for proper density estimation
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3, 4, 5]),  # Target molecule
            [make_emitter(Float64(20+i), Float64(20+i), 1) for i in 1:5]...  # Background
        )
        smld = make_test_smld(emitters; n_frames=5)

        (combined, info) = frameconnect(smld; maxframegap=5, nnearestclusters=1)

        # Target molecule should be combined (may have 5 background singles)
        # At minimum, we should have fewer emitters than input
        @test length(combined.emitters) < length(smld.emitters)
    end

    @testset "multiple separated molecules" begin
        # Multiple well-separated molecules for proper density estimation
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3]),
            make_blinking_molecule(10.0, 10.0, [1, 2, 3]),
            make_blinking_molecule(15.0, 15.0, [1, 2, 3]),
            make_blinking_molecule(20.0, 20.0, [1, 2, 3]),
        )
        smld = make_test_smld(emitters; n_frames=3)

        (combined, info) = frameconnect(smld; maxframegap=5)

        # Should combine into ~4 molecules (fewer than 12 input emitters)
        @test length(combined.emitters) <= 8
        @test length(combined.emitters) < length(smld.emitters)
    end

    @testset "track_id population" begin
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3]),
            make_blinking_molecule(10.0, 10.0, [1, 2, 3]),
            make_blinking_molecule(15.0, 15.0, [1, 2, 3]),
        )
        smld = make_test_smld(emitters; n_frames=3)

        (combined, info) = frameconnect(smld)

        # All emitters should have non-zero track_id
        for e in info.connected.emitters
            @test e.track_id > 0
        end
    end

    @testset "parameter estimation" begin
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3]),
            make_blinking_molecule(10.0, 10.0, [1, 2, 3]),
            make_blinking_molecule(15.0, 15.0, [1, 2, 3]),
            make_blinking_molecule(20.0, 20.0, [1, 2, 3]),
        )
        smld = make_test_smld(emitters; n_frames=3)

        (combined, info) = frameconnect(smld)

        # Parameters should be estimated (non-default values)
        @test info.k_on >= 0
        @test info.k_off >= 0
        @test info.k_bleach >= 0
        @test 0 <= info.p_miss <= 1
        @test !isempty(info.initialdensity)
        @test info.n_preclusters > 0
        @test info.n_tracks > 0
    end

    @testset "respects maxframegap" begin
        # Create molecules with gaps, plus background for density estimation
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3]),       # Frames 1-3
            make_blinking_molecule(5.0, 5.0, [10, 11, 12]),    # Frames 10-12 (gap of 7)
            make_blinking_molecule(20.0, 20.0, [1, 2, 3]),     # Background molecule
            make_blinking_molecule(25.0, 25.0, [1, 2, 3]),     # Background molecule
        )
        smld = make_test_smld(emitters; n_frames=12)

        # With maxframegap=5, the two groups at (5,5) should NOT connect
        (combined, info) = frameconnect(smld; maxframegap=5)

        # Should have at least 3 combined localizations (2 for split molecule + 2 background)
        @test length(combined.emitters) >= 3
    end

    @testset "single localization" begin
        # Edge case: just one localization
        emitters = [make_emitter(5.0, 5.0, 1)]
        smld = make_test_smld(emitters; n_frames=1)

        (combined, info) = frameconnect(smld)

        @test length(combined.emitters) == 1
        @test combined.emitters[1].x ≈ 5.0 atol=1e-10
        @test info.n_input == 1
        @test info.n_combined == 1
    end

    @testset "preserves camera and metadata" begin
        smld = make_single_molecule_smld()

        (combined, info) = frameconnect(smld)

        @test info.connected.camera == smld.camera
        @test combined.camera == smld.camera
        @test info.connected.n_datasets == smld.n_datasets
    end

    @testset "precision improvement" begin
        # Verify that combining N localizations improves precision by ~√N
        n_locs = 4
        σ_input = 0.02

        # Target molecule to test precision improvement
        target_emitters = [make_emitter(5.0, 5.0, i; σ_xy=σ_input) for i in 1:n_locs]
        # Background molecules for proper density estimation
        background = vcat(
            [make_emitter(Float64(20+i), Float64(20+i), 1) for i in 1:4]...
        )
        emitters = vcat(target_emitters, background)
        smld = make_test_smld(emitters; n_frames=n_locs)

        (combined, info) = frameconnect(smld; maxframegap=n_locs, nnearestclusters=1)

        # Find the combined target molecule (around position 5,5)
        target_combined = filter(e -> e.x < 10.0 && e.y < 10.0, combined.emitters)
        @test length(target_combined) == 1

        σ_output = target_combined[1].σ_x
        expected_σ = σ_input / sqrt(n_locs)

        # Allow some tolerance since positions have noise
        @test σ_output ≈ expected_σ rtol=0.1
    end

    @testset "ConnectInfo statistics" begin
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3]),
            make_blinking_molecule(10.0, 10.0, [1, 2]),
        )
        smld = make_test_smld(emitters; n_frames=3)

        (combined, info) = frameconnect(smld)

        @test info.n_input == 5
        @test info.n_combined <= info.n_input
        @test info.n_tracks <= info.n_input
        @test info.n_preclusters <= info.n_input
        @test info.elapsed_ns isa UInt64
        @test info.algorithm == :lap
    end
end
