using Test
using SMLMFrameConnection
using SMLMData

# Load test fixtures
include("test_helpers.jl")

@testset "SMLMFrameConnection.jl" begin
    include("test_types.jl")
    include("test_combinelocalizations.jl")
    include("test_frameconnect.jl")
    include("test_defineidealFC.jl")
    include("test_calibration.jl")
end
