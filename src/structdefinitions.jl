# This file defines some struct types used in the FrameConnection package.

"""
    ParamStruct(initialdensity::Vector{Float32}, 
        k_on::Float32, k_off::Float32, k_bleach::Float32, p_miss::Float32, 
        nsigmadev::Float32, maxframegap::Int, nnearestclusters::Int)

Structure of parameters needed for frame-connection.

# Fields
-`initialdensity`: Density of emitters at the start of the experiment.
                   (see estimatedensities()) (emitters/pixel^2)
-`nnearestclusters`: Number of nearest preclusters used for local density 
                     estimates. (default = 2)(see estimatedensities())
-`k_on`: Rate at which dark emitters convert to the visible state. 
         (see estimateparams())(1/frame)
-`k_off`: Rate at which visible emitters are converted to the reversible dark
          state. (see estimateparams())(1/frame)
-`k_bleach`: Rate at which visible emitters are irreversibly photobleached.
             (see estimateparams())(1/frame)
-`p_miss`: Probability of missing a localization of a visible emitter.
-`nsigmadev`: Multiplier of localization errors that defines a pre-clustering
              distance threshold. (default = 5)(see precluster())(pixels)
-`maxframegap`: Maximum frame gap between temporally adjacent localizations in
                a precluster. (default = 5)(see precluster())(frames)
-`nmaxnn`: Maximum number of nearest-neighbors inspected for precluster 
           membership.  Ideally, this would be set to inf, but that's not
           feasible for most data. (default = 2)(see precluster())
        

"""
mutable struct ParamStruct
    initialdensity::Vector{Float32}
    nnearestclusters::Int
    k_on::Float32
    k_off::Float32
    k_bleach::Float32
    p_miss::Float32
    nsigmadev::Float32
    maxframegap::Int
    nmaxnn::Int
end
ParamStruct() = ParamStruct([], 2, 0, 0, 0, 0, 5, 5, 2)