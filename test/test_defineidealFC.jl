@testset "defineidealFC" begin
    @testset "return types" begin
        smld = make_prelabeled_smld()
        smld_connected, smld_combined = defineidealFC(smld)

        @test smld_connected isa BasicSMLD
        @test smld_combined isa BasicSMLD
    end

    @testset "groups by track_id" begin
        # Create data where track_id represents true emitter identity
        emitters = [
            make_emitter(5.0, 5.0, 1; track_id=1),
            make_emitter(5.01, 5.01, 2; track_id=1),
            make_emitter(5.02, 4.99, 3; track_id=1),
            make_emitter(10.0, 10.0, 1; track_id=2),
            make_emitter(10.01, 10.01, 2; track_id=2),
        ]
        smld = make_test_smld(emitters; n_frames=3)

        smld_connected, smld_combined = defineidealFC(smld; max_frame_gap=5)

        # Should have 2 combined localizations (one per unique emitter)
        @test length(smld_combined.emitters) == 2
    end

    @testset "respects max_frame_gap" begin
        # Same emitter (track_id=1) but with frame gap > max_frame_gap
        emitters = [
            make_emitter(5.0, 5.0, 1; track_id=1),
            make_emitter(5.0, 5.0, 2; track_id=1),
            make_emitter(5.0, 5.0, 10; track_id=1),  # Gap of 8 frames
            make_emitter(5.0, 5.0, 11; track_id=1),
        ]
        smld = make_test_smld(emitters; n_frames=11)

        # With max_frame_gap=5, should split into 2 blinking events
        smld_connected, smld_combined = defineidealFC(smld; max_frame_gap=5)

        @test length(smld_combined.emitters) == 2
    end

    @testset "consecutive frames connect" begin
        # All consecutive frames should connect
        emitters = [make_emitter(5.0, 5.0, i; track_id=1) for i in 1:5]
        smld = make_test_smld(emitters; n_frames=5)

        _, smld_combined = defineidealFC(smld; max_frame_gap=5)

        @test length(smld_combined.emitters) == 1
    end

    @testset "different emitters stay separate" begin
        # Two different emitters (different track_id) at same location
        emitters = [
            make_emitter(5.0, 5.0, 1; track_id=1),
            make_emitter(5.0, 5.0, 2; track_id=2),
            make_emitter(5.0, 5.0, 3; track_id=1),
            make_emitter(5.0, 5.0, 4; track_id=2),
        ]
        smld = make_test_smld(emitters; n_frames=4)

        _, smld_combined = defineidealFC(smld; max_frame_gap=5)

        # Should have 2 combined localizations (one per emitter identity)
        @test length(smld_combined.emitters) == 2
    end

    @testset "single localization" begin
        emitters = [make_emitter(5.0, 5.0, 1; track_id=1)]
        smld = make_test_smld(emitters; n_frames=1)

        smld_connected, smld_combined = defineidealFC(smld)

        @test length(smld_combined.emitters) == 1
        @test smld_combined.emitters[1].x ≈ 5.0 atol=1e-10
    end

    @testset "connected track_id compressed" begin
        # After connection, track_ids should be 1:N
        emitters = [
            make_emitter(5.0, 5.0, 1; track_id=100),
            make_emitter(5.0, 5.0, 2; track_id=100),
            make_emitter(10.0, 10.0, 1; track_id=200),
        ]
        smld = make_test_smld(emitters; n_frames=2)

        smld_connected, _ = defineidealFC(smld)

        track_ids = sort(unique(e.track_id for e in smld_connected.emitters))
        # Should be compressed to 1, 2 (not 100, 200)
        @test track_ids == [1, 2]
    end
end
