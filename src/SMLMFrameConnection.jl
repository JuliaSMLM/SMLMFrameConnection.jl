module SMLMFrameConnection

using SMLMData
using SMLMData: BasicSMLD, Emitter2DFit
using Hungarian
using NearestNeighbors
using Optim
using StatsBase
using Statistics
using StructArrays


include("structdefinitions.jl")
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