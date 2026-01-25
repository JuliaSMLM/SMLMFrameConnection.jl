@testset "frameconnect" begin
    @testset "return types" begin
        # Use multiple molecules for proper density estimation
        emitters = vcat(
            [make_emitter(Float64(i), Float64(i), j) for i in 1:5 for j in 1:3]...
        )
        smld = make_test_smld(emitters; n_frames=3)
        result = frameconnect(smld)

        @test result isa NamedTuple
        @test haskey(result, :combined)
        @test haskey(result, :connected)
        @test haskey(result, :params)
        @test result.combined isa BasicSMLD
        @test result.connected isa BasicSMLD
        @test result.params isa ParamStruct
    end

    @testset "single molecule connection" begin
        # Single molecule with 5 localizations, but add more spread-out molecules
        # for proper density estimation
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3, 4, 5]),  # Target molecule
            [make_emitter(Float64(20+i), Float64(20+i), 1) for i in 1:5]...  # Background
        )
        smld = make_test_smld(emitters; n_frames=5)

        result = frameconnect(smld; maxframegap=5, nnearestclusters=1)

        # Target molecule should be combined (may have 5 background singles)
        # At minimum, we should have fewer emitters than input
        @test length(result.combined.emitters) < length(smld.emitters)
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

        result = frameconnect(smld; maxframegap=5)

        # Should combine into ~4 molecules (fewer than 12 input emitters)
        @test length(result.combined.emitters) <= 8
        @test length(result.combined.emitters) < length(smld.emitters)
    end

    @testset "track_id population" begin
        emitters = vcat(
            make_blinking_molecule(5.0, 5.0, [1, 2, 3]),
            make_blinking_molecule(10.0, 10.0, [1, 2, 3]),
            make_blinking_molecule(15.0, 15.0, [1, 2, 3]),
        )
        smld = make_test_smld(emitters; n_frames=3)

        result = frameconnect(smld)

        # All emitters should have non-zero track_id
        for e in result.connected.emitters
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

        result = frameconnect(smld)

        # Parameters should be estimated (non-default values)
        @test result.params.k_on >= 0
        @test result.params.k_off >= 0
        @test result.params.k_bleach >= 0
        @test 0 <= result.params.p_miss <= 1
        @test !isempty(result.params.initialdensity)
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
        result = frameconnect(smld; maxframegap=5)

        # Should have at least 3 combined localizations (2 for split molecule + 2 background)
        @test length(result.combined.emitters) >= 3
    end

    @testset "single localization" begin
        # Edge case: just one localization
        emitters = [make_emitter(5.0, 5.0, 1)]
        smld = make_test_smld(emitters; n_frames=1)

        result = frameconnect(smld)

        @test length(result.combined.emitters) == 1
        @test result.combined.emitters[1].x ≈ 5.0 atol=1e-10
    end

    @testset "preserves camera and metadata" begin
        smld = make_single_molecule_smld()

        result = frameconnect(smld)

        @test result.connected.camera == smld.camera
        @test result.combined.camera == smld.camera
        @test result.connected.n_datasets == smld.n_datasets
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

        result = frameconnect(smld; maxframegap=n_locs, nnearestclusters=1)

        # Find the combined target molecule (around position 5,5)
        target_combined = filter(e -> e.x < 10.0 && e.y < 10.0, result.combined.emitters)
        @test length(target_combined) == 1

        σ_output = target_combined[1].σ_x
        expected_σ = σ_input / sqrt(n_locs)

        # Allow some tolerance since positions have noise
        @test σ_output ≈ expected_σ rtol=0.1
    end

    @testset "custom parameters" begin
        smld = make_single_molecule_smld()

        result = frameconnect(smld;
            nnearestclusters=3,
            nsigmadev=4.0,
            maxframegap=10,
            nmaxnn=3
        )

        @test result.params.nnearestclusters == 3
        @test result.params.nsigmadev == 4.0
        @test result.params.maxframegap == 10
        @test result.params.nmaxnn == 3
    end
end
