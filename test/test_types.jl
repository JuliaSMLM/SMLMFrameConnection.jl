@testset "Types" begin
    @testset "ParamStruct" begin
        @testset "default constructor" begin
            params = ParamStruct()
            @test params isa ParamStruct
            @test params.initialdensity == []
            @test params.nnearestclusters == 2
            @test params.k_on == 0.0
            @test params.k_off == 0.0
            @test params.k_bleach == 0.0
            @test params.p_miss == 0.0
            @test params.nsigmadev == 5.0
            @test params.maxframegap == 5
            @test params.nmaxnn == 2
        end

        @testset "field mutation" begin
            params = ParamStruct()

            # Should be mutable
            params.k_on = 0.1
            @test params.k_on == 0.1

            params.k_off = 0.5
            @test params.k_off == 0.5

            params.initialdensity = [1.0, 2.0]
            @test params.initialdensity == [1.0, 2.0]
        end

        @testset "full constructor" begin
            params = ParamStruct([1.0], 3, 0.1, 0.5, 0.01, 0.05, 4.0, 10, 3)
            @test params.initialdensity == [1.0]
            @test params.nnearestclusters == 3
            @test params.k_on == 0.1
            @test params.k_off == 0.5
            @test params.k_bleach == 0.01
            @test params.p_miss == 0.05
            @test params.nsigmadev == 4.0
            @test params.maxframegap == 10
            @test params.nmaxnn == 3
        end
    end
end
