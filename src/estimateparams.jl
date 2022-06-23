using SMLMData
using Optim
using Statistics

"""
    k_on, k_off, k_bleach, p_miss, nemitters = 
        estimateparams(smld::SMLMData.SMLD2D, clusterdata::Vector{Matrix{Float32}})

Estimate rate parameters from the clusters in `smld` and `clusterdata`.

# Description
The rate parameters `k_on`, `k_off`, and `k_bleach`, the probability of missing
a localization for a visible emitter `p_miss`, and the underlying number of
emitters `nemitters` are estimated from the pre-clusters in 
`smld`/`clusterdata`.  This is done by assuming that pre-clusters are, on 
average, representative of a single blinking event of a single emitter.

# Citation
Schodt David J., Lidke Keith A., "Spatiotemporal Clustering of Repeated 
Super-Resolution Localizations via Linear Assignment Problem", 
Front. Bioinform., 20 October 2021
https://doi.org/10.3389/fbinf.2021.724325
"""
function estimateparams(smld::SMLMData.SMLD2D, 
    clusterdata::Vector{Matrix{Float32}})

    # Compute some quantities needed from each of the preclusters.
    nclusters = size(clusterdata)[1]
    clusterdurations, nobservations = computeclusterinfo(clusterdata)

    # Estimate the sum of the off rate and the bleaching rate, as well as the
    # probability of missing a localization.
    k_offpbleach = -log(1 - 1 / Statistics.mean(clusterdurations))
    k_offpbleach = isinf(k_offpbleach) ? 1.0 : k_offpbleach
    p_miss = 1 - Statistics.mean(nobservations ./ clusterdurations)

    # Compute other parameters by fitting the observed number of localizations
    # over time (see citation for description of this math).
    nlocperframe = counts(smld.framenum)
    nlocperframe = nlocperframe[nlocperframe.!=0]
    frames = convert(Array{Float32}, unique(smld.framenum))
    nloccumulative = cumsum(nlocperframe)
    k_off = nclusters / nloccumulative[end]
    k_bleach = maximum([10^-5 k_offpbleach - k_off])
    k(k_on) = k_on + k_offpbleach
    l1(k_on) = k_on * k_bleach / k(k_on)
    l2(k_on) = k(k_on) - l1(k_on)
    model(x) = ceil(x[1]) * (1.0 .- p_miss) * (x[2] / k(x[2])) *
               ((1 / l1(x[2])) * (1.0 .- exp.(-l1(x[2]) * (frames .- 1.0))) .-
                (1 / l2(x[2])) * (1.0 .- exp.(-l2(x[2]) * (frames .- 1.0))))
    costfunction(x) = Statistics.mean((nloccumulative .- model(x)) .^ 2)
    lowerbound = [maximum(nlocperframe) 10^-5]
    upperbound = [nclusters nloccumulative[end] / frames[end]]
    x1_guess = ceil(nclusters * k_bleach / (k_off*(1.0-p_miss))) # DJS 22/06/23: better guess than suggested in paper
    x1_guess = (x1_guess>lowerbound[1]) & (x1_guess<upperbound[1]) ?
               x1_guess : (upperbound[1]-lowerbound[1]) / 2.0
    initialguess = [x1_guess 1.0 / frames[end]]
    inner_optimizer = GradientDescent()
    optimresults = optimize(costfunction, lowerbound, upperbound, initialguess,
        Fminbox(inner_optimizer))
    xhat = optimresults.minimizer
    nemitters = xhat[1]
    k_on = xhat[2]

    return k_on, k_off, k_bleach, p_miss, nemitters
end