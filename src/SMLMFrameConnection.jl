module SMLMFrameConnection

using SMLMData
using Hungarian
using NearestNeighbors
using Optim
using StatsBase
using Statistics

export frameconnect, defineidealFC, combinelocalizations
export ConnectConfig, ConnectInfo, ParamStruct

include("structdefinitions.jl")
include("connectinfo.jl")
include("precluster.jl")
include("defineidealFC.jl")
include("organizeclusters.jl")
include("computeclusterinfo.jl")
include("estimateparams.jl")
include("estimatedensities.jl")
include("create_costmatrix.jl")
include("solveLAP.jl")
include("linkclusters.jl")
include("compress_connectID.jl")
include("combinelocalizations.jl")
include("connectlocalizations.jl")
include("frameconnect.jl")


end