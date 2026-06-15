# Robustness tests for empty (0) and tiny (1-2) localization inputs.
#
# frameconnect and its building blocks previously crashed on near-empty data
# (e.g. `MethodError: zero(Type{Any})` / `reducing over an empty collection`)
# when reductions like `counts`/`mean`/`extrema` ran over empty collections.
# These tests pin the graceful behavior: a passthrough/empty result, no crash.

@testset "empty and tiny inputs" begin
    cam = make_test_camera()
    # Build SMLDs directly here: make_test_smld reduces over emitters (maximum
    # frame) and so cannot construct a 0-emitter SMLD.
    empty_smld(; nframes=0) =
        BasicSMLD(Emitter2DFit{Float64}[], cam, nframes, 1, Dict{String,Any}())

    @testset "frameconnect: empty SMLD (n_frames=0)" begin
        smld = empty_smld()
        (combined, info) = frameconnect(smld)

        @test combined isa BasicSMLD
        @test isempty(combined.emitters)
        @test info isa FrameConnectInfo{Float64}
        @test info.n_input == 0
        @test info.n_tracks == 0
        @test info.n_combined == 0
        @test info.n_preclusters == 0
        @test info.n_filtered == 0
        @test info.calibration === nothing
        @test info.elapsed_s >= 0
    end

    @testset "frameconnect: empty SMLD (n_frames>0)" begin
        smld = empty_smld(nframes=100)
        (combined, info) = frameconnect(smld)
        @test isempty(combined.emitters)
        @test info.n_combined == 0
    end

    @testset "frameconnect: empty SMLD preserves camera/datasets" begin
        smld = empty_smld()
        (combined, info) = frameconnect(smld)
        @test combined.camera == smld.camera
        @test info.connected.camera == smld.camera
        @test combined.n_datasets == smld.n_datasets
    end

    @testset "frameconnect: empty SMLD with calibration falls back gracefully" begin
        # Calibration enabled + empty input must mirror the tiny-input contract:
        # a CalibrationResult with calibration_applied=false, NOT `nothing`, so
        # downstream code reading info.calibration.* sees a uniform shape.
        smld = empty_smld()
        config = FrameConnectConfig(calibration=CalibrationConfig())
        (combined, info) = frameconnect(smld, config)

        @test isempty(combined.emitters)
        @test info.calibration !== nothing
        @test info.calibration.calibration_applied == false
        @test info.calibration.n_pairs == 0
    end

    @testset "frameconnect: single localization" begin
        smld = make_test_smld([make_emitter(5.0, 5.0, 1)]; n_frames=1)
        (combined, info) = frameconnect(smld)
        @test length(combined.emitters) == 1
        @test info.n_input == 1
        @test info.n_combined == 1
    end

    @testset "frameconnect: single localization with calibration" begin
        smld = make_test_smld([make_emitter(5.0, 5.0, 1)]; n_frames=1)
        config = FrameConnectConfig(calibration=CalibrationConfig())
        (combined, info) = frameconnect(smld, config)
        @test length(combined.emitters) == 1
        @test info.calibration !== nothing
        @test info.calibration.calibration_applied == false  # 0 frame-to-frame pairs
    end

    @testset "frameconnect: two localizations (one molecule)" begin
        # Small jitter (via make_blinking_molecule) avoids a degenerate
        # zero-area cluster from two identical positions.
        emitters = make_blinking_molecule(5.0, 5.0, [1, 2])
        smld = make_test_smld(emitters; n_frames=2)
        (combined, info) = frameconnect(smld)
        @test combined isa BasicSMLD
        @test info.n_input == 2
        @test length(combined.emitters) >= 1
        @test length(combined.emitters) <= 2
    end

    @testset "frameconnect: two localizations (separated)" begin
        emitters = [make_emitter(5.0, 5.0, 1), make_emitter(8.0, 8.0, 1)]
        smld = make_test_smld(emitters; n_frames=1)
        (combined, info) = frameconnect(smld)
        @test info.n_input == 2
        @test length(combined.emitters) >= 1
    end

    @testset "frameconnect: empty SMLD (abstract emitter eltype)" begin
        # An untyped/abstract emitter vector makes the internal reductions
        # produce Any[] collections, which is what triggered the original
        # `MethodError: zero(Type{Any})`. The empty guard must handle this
        # regardless of element type.
        smld = BasicSMLD(Emitter2DFit[], cam, 0, 1, Dict{String,Any}())
        (combined, info) = frameconnect(smld)
        @test isempty(combined.emitters)
        @test info.n_input == 0
        @test info.n_combined == 0
    end

    @testset "combinelocalizations: empty passthrough" begin
        smld = empty_smld()
        out = combinelocalizations(smld)
        @test out isa BasicSMLD
        @test isempty(out.emitters)
    end

    @testset "precluster: empty passthrough" begin
        smld = empty_smld()
        out = SMLMFrameConnection.precluster(smld)
        @test out isa BasicSMLD
        @test isempty(out.emitters)
    end
end
