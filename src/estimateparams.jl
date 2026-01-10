using SMLMData
using Optim
using Statistics

"""
    k_on, k_off, k_bleach, p_miss, nemitters =
        estimateparams(smld::BasicSMLD{T,SMLMData.Emitter2DFit{T}},
                       clusterdata::Vector{Matrix{Float32}}) where T

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
function estimateparams(smld::BasicSMLD{T,SMLMData.Emitter2DFit{T}},
    clusterdata::Vector{Matrix{Float32}}) where T

    # Extract frame numbers from emitters
    framenum = [e.frame for e in smld.emitters]

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
    nlocperframe = counts(framenum)
    nlocperframe = nlocperframe[nlocperframe.!=0]
    frames = convert(Array{Float64}, unique(framenum))
    nloccumulative = cumsum(nlocperframe)
    k_off = nclusters / nloccumulative[end]
    k_bleach = maximum([10^-5 k_offpbleach - k_off])

    # Pre-compute frames-1 once (avoid repeated allocation)
    frames_m1 = frames .- 1.0
    nframes = length(frames)

    # Allocation-free cost function for optimization
    # Uses a loop instead of broadcasting to avoid temporary arrays
    function costfunction(x)
        nem, kon = x[1], x[2]
        k_val = kon + k_offpbleach
        l1 = kon * k_bleach / k_val
        l2 = k_val - l1
        coef = ceil(nem) * (1.0 - p_miss) * (kon / k_val)

        # Guard against division by zero
        if l1 ≈ 0.0 || l2 ≈ 0.0
            return 1e10
        end

        inv_l1 = 1.0 / l1
        inv_l2 = 1.0 / l2

        mse = 0.0
        @inbounds for i in 1:nframes
            fm1 = frames_m1[i]
            pred = coef * (inv_l1 * (1.0 - exp(-l1 * fm1)) -
                          inv_l2 * (1.0 - exp(-l2 * fm1)))
            diff = nloccumulative[i] - pred
            mse += diff * diff
        end
        return mse / nframes
    end

    lowerbound = [Float64(maximum(nlocperframe)); 1e-5]
    upperbound = [Float64(nclusters); nloccumulative[end] / frames[end]]

    # Ensure upperbound > lowerbound
    if upperbound[1] <= lowerbound[1]
        upperbound[1] = lowerbound[1] + 1.0
    end
    if upperbound[2] <= lowerbound[2]
        upperbound[2] = lowerbound[2] + 1e-4
    end

    x1_guess = ceil(nclusters * k_bleach / (k_off*(1.0-p_miss))) # DJS 22/06/23: better guess than suggested in paper
    x1_guess = (x1_guess>lowerbound[1]) && (x1_guess<upperbound[1]) ?
               x1_guess : (upperbound[1]+lowerbound[1]) / 2.0
    initialguess = [x1_guess; 1.0 / frames[end]]

    # Ensure initialguess is within bounds
    initialguess[1] = clamp(initialguess[1], lowerbound[1], upperbound[1])
    initialguess[2] = clamp(initialguess[2], lowerbound[2], upperbound[2])

    # Use NelderMead: derivative-free, fast convergence for 2-parameter problems
    # (Previous: GradientDescent was extremely slow due to linear convergence
    # and finite-difference gradient computation)
    optimresults = optimize(costfunction, lowerbound, upperbound, initialguess,
        Fminbox(NelderMead()), Optim.Options(iterations=500))
    xhat = optimresults.minimizer
    nemitters = xhat[1]
    k_on = xhat[2]

    return k_on, k_off, k_bleach, p_miss, nemitters
end
