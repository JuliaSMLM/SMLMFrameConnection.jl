@testset "combinelocalizations" begin
    @testset "single localization per cluster" begin
        # Each emitter has unique track_id - nothing to combine
        emitters = [
            make_emitter(5.0, 5.0, 1; track_id=1),
            make_emitter(10.0, 10.0, 2; track_id=2),
            make_emitter(15.0, 15.0, 3; track_id=3),
        ]
        smld = make_test_smld(emitters)
        result = combinelocalizations(smld)

        @test length(result.emitters) == 3
        # Positions should be unchanged
        @test result.emitters[1].x ≈ 5.0 atol=1e-10
        @test result.emitters[2].x ≈ 10.0 atol=1e-10
        @test result.emitters[3].x ≈ 15.0 atol=1e-10
    end

    @testset "two localizations same track_id" begin
        # Two localizations at same position with same uncertainty
        σ = 0.02
        emitters = [
            make_emitter(5.0, 5.0, 1; σ_pos=σ, track_id=1, photons=1000.0),
            make_emitter(5.0, 5.0, 2; σ_pos=σ, track_id=1, photons=1200.0),
        ]
        smld = make_test_smld(emitters)
        result = combinelocalizations(smld)

        @test length(result.emitters) == 1

        e = result.emitters[1]
        # Position should be at 5.0 (both inputs identical)
        @test e.x ≈ 5.0 atol=1e-10
        @test e.y ≈ 5.0 atol=1e-10

        # Combined uncertainty: σ_combined = σ / √2
        expected_σ = σ / sqrt(2)
        @test e.σ_x ≈ expected_σ atol=1e-10
        @test e.σ_y ≈ expected_σ atol=1e-10

        # Photons should be summed
        @test e.photons ≈ 2200.0 atol=1e-10

        # track_id preserved
        @test e.track_id == 1
    end

    @testset "MLE weighted mean" begin
        # Two localizations with different uncertainties
        # Weighted mean should favor lower uncertainty
        emitters = [
            make_emitter(5.0, 5.0, 1; σ_pos=0.01, track_id=1),  # High precision
            make_emitter(5.1, 5.1, 2; σ_pos=0.10, track_id=1),  # Low precision
        ]
        smld = make_test_smld(emitters)
        result = combinelocalizations(smld)

        @test length(result.emitters) == 1

        e = result.emitters[1]
        # Result should be closer to (5.0, 5.0) due to higher weight
        # Weight ratio: (1/0.01²) : (1/0.1²) = 10000 : 100 = 100:1
        # Weighted mean: (5.0*100 + 5.1*1) / 101 ≈ 5.001
        @test e.x ≈ 5.0 atol=0.01
        @test e.y ≈ 5.0 atol=0.01
    end

    @testset "three localizations" begin
        σ = 0.02
        emitters = [
            make_emitter(5.00, 5.00, 1; σ_pos=σ, track_id=1),
            make_emitter(5.02, 5.02, 2; σ_pos=σ, track_id=1),
            make_emitter(4.98, 4.98, 3; σ_pos=σ, track_id=1),
        ]
        smld = make_test_smld(emitters)
        result = combinelocalizations(smld)

        @test length(result.emitters) == 1

        e = result.emitters[1]
        # Mean position (equal weights)
        @test e.x ≈ 5.0 atol=0.01
        @test e.y ≈ 5.0 atol=0.01

        # Combined uncertainty: σ / √3
        expected_σ = σ / sqrt(3)
        @test e.σ_x ≈ expected_σ atol=1e-10
        @test e.σ_y ≈ expected_σ atol=1e-10
    end

    @testset "multiple clusters" begin
        smld = make_prelabeled_smld()
        result = combinelocalizations(smld)

        # Should have 3 clusters (track_id 1, 2, 3)
        @test length(result.emitters) == 3

        # Sort by track_id to check each cluster
        sorted = sort(result.emitters, by=e->e.track_id)

        # Cluster 1: 2 localizations
        @test sorted[1].track_id == 1

        # Cluster 2: 3 localizations
        @test sorted[2].track_id == 2

        # Cluster 3: 1 localization (unchanged)
        @test sorted[3].track_id == 3
        @test sorted[3].x ≈ 15.0 atol=1e-10
    end

    @testset "photon summation" begin
        emitters = [
            make_emitter(5.0, 5.0, 1; track_id=1, photons=1000.0),
            make_emitter(5.0, 5.0, 2; track_id=1, photons=1500.0),
            make_emitter(5.0, 5.0, 3; track_id=1, photons=2000.0),
        ]
        smld = make_test_smld(emitters)
        result = combinelocalizations(smld)

        @test length(result.emitters) == 1
        @test result.emitters[1].photons ≈ 4500.0 atol=1e-10
    end

    @testset "Float32 precision at SMLM scale" begin
        # Realistic SMLM values: σ ~ 6nm = 0.006μm in Float32
        # This tests that precision-weighted combination doesn't lose
        # accuracy due to Float32 determinant underflow
        σ = Float32(0.006)  # 6 nm
        x0 = Float32(25.0)
        y0 = Float32(25.0)

        emitters = [
            Emitter2DFit{Float32}(
                x0 + Float32(0.001), y0 - Float32(0.001), 1000f0, 10f0,
                σ, σ, 0f0, 31.6f0, 3.16f0, 1, 1, 1, 0),
            Emitter2DFit{Float32}(
                x0 - Float32(0.001), y0 + Float32(0.001), 1200f0, 12f0,
                σ, σ, 0f0, 34.6f0, 3.46f0, 2, 1, 1, 0),
            Emitter2DFit{Float32}(
                x0, y0, 1100f0, 11f0,
                σ, σ, 0f0, 33.2f0, 3.32f0, 3, 1, 1, 0),
        ]

        smld = make_test_smld(emitters; n_frames=3)
        result = combinelocalizations(smld)

        @test length(result.emitters) == 1
        e = result.emitters[1]

        # Combined σ must be LESS than individual σ (σ/√3 for equal weights)
        @test e.σ_x < σ
        @test e.σ_y < σ

        # Expected: σ/√3 ≈ 0.00346
        expected_σ = σ / sqrt(Float32(3))
        @test e.σ_x ≈ expected_σ atol=1e-5
        @test e.σ_y ≈ expected_σ atol=1e-5

        # Position near centroid
        @test e.x ≈ x0 atol=0.002
        @test e.y ≈ y0 atol=0.002
    end

    @testset "preserves metadata" begin
        smld = make_prelabeled_smld()
        result = combinelocalizations(smld)

        @test result.camera == smld.camera
        @test result.n_frames == smld.n_frames
        @test result.n_datasets == smld.n_datasets
    end
end
