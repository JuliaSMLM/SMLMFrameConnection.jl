@testset "Types" begin
    @testset "FrameConnectInfo" begin
        # Create minimal test data
        emitters = [SMLMData.Emitter2DFit{Float64}(
            5.0, 5.0, 1000.0, 10.0, 0.02, 0.02, 0.0, 10.0, 1.0, 1, 1, 1, 1
        )]
        camera = SMLMData.IdealCamera(1:64, 1:64, 0.1)
        smld = BasicSMLD(emitters, camera, 1, 1, Dict{String,Any}())

        info = FrameConnectInfo{Float64}(
            smld, 10, 5, 5, 0.1, 0.5, 0.01, 0.05, [1.0, 2.0], 1.5, :lap, 3, nothing
        )

        @test info isa FrameConnectInfo{Float64}
        @test info.connected === smld
        @test info.n_input == 10
        @test info.n_tracks == 5
        @test info.n_combined == 5
        @test info.k_on == 0.1
        @test info.k_off == 0.5
        @test info.k_bleach == 0.01
        @test info.p_miss == 0.05
        @test info.initial_density == [1.0, 2.0]
        @test info.elapsed_s == 1.5
        @test info.algorithm == :lap
        @test info.n_preclusters == 3
        @test info.calibration === nothing
    end

    @testset "ParamStruct" begin
        @testset "default constructor" begin
            params = ParamStruct()
            @test params isa ParamStruct
            @test params.initial_density == []
            @test params.n_density_neighbors == 2
            @test params.k_on == 0.0
            @test params.k_off == 0.0
            @test params.k_bleach == 0.0
            @test params.p_miss == 0.0
            @test params.max_sigma_dist == 5.0
            @test params.max_frame_gap == 5
            @test params.max_neighbors == 2
        end

        @testset "field mutation" begin
            params = ParamStruct()

            # Should be mutable
            params.k_on = 0.1
            @test params.k_on == 0.1

            params.k_off = 0.5
            @test params.k_off == 0.5

            params.initial_density = [1.0, 2.0]
            @test params.initial_density == [1.0, 2.0]
        end

        @testset "full constructor" begin
            params = ParamStruct([1.0], 3, 0.1, 0.5, 0.01, 0.05, 4.0, 10, 3)
            @test params.initial_density == [1.0]
            @test params.n_density_neighbors == 3
            @test params.k_on == 0.1
            @test params.k_off == 0.5
            @test params.k_bleach == 0.01
            @test params.p_miss == 0.05
            @test params.max_sigma_dist == 4.0
            @test params.max_frame_gap == 10
            @test params.max_neighbors == 3
        end
    end
end
