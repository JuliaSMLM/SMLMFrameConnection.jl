using Optim
using Statistics
using StatsBase

"""
    estimateparams(smld::BasicSMLD, clusterdata::Vector{Matrix{Float32}})

Estimate rate parameters from the clusters in `smld` and `clusterdata`.

# Description
The rate parameters `k_on`, `k_off`, and `k_bleach`, the probability of missing
a localization for a visible emitter `p_miss`, and the underlying number of
emitters `nemitters` are estimated from the pre-clusters in `smld`/`clusterdata`.
This is done by assuming that pre-clusters are, on average, representative of a
single blinking event of a single emitter.

# Returns
Tuple of (k_on, k_off, k_bleach, p_miss, nemitters)

# Citation
Schodt David J., Lidke Keith A., "Spatiotemporal Clustering of Repeated
Super-Resolution Localizations via Linear Assignment Problem",
Front. Bioinform., 20 October 2021
https://doi.org/10.3389/fbinf.2021.724325
"""
function estimateparams(
    smld::BasicSMLD{T, E},
    clusterdata::Vector{Matrix{Float32}}
) where {T, E<:Emitter2DFit}
    # Compute statistics from preclusters
    nclusters = length(clusterdata)
    clusterdurations, nobservations = computeclusterinfo(clusterdata)

    # Estimate sum of off rate and bleaching rate, and miss probability
    k_offpbleach = -log(1 - 1 / mean(clusterdurations))
    k_offpbleach = isinf(k_offpbleach) ? 1.0 : k_offpbleach
    p_miss = 1 - mean(nobservations ./ clusterdurations)

    # Extract frame numbers from emitters
    frames_all = [e.frame for e in smld.emitters]
    nlocperframe = counts(frames_all)
    nlocperframe = nlocperframe[nlocperframe .!= 0]
    frames = Float32.(sort(unique(frames_all)))
    nloccumulative = cumsum(nlocperframe)

    # Initial parameter estimates
    k_off = nclusters / nloccumulative[end]
    k_bleach = max(1e-5, k_offpbleach - k_off)

    # Helper functions for fitting model
    k(k_on) = k_on + k_offpbleach
    l1(k_on) = k_on * k_bleach / k(k_on)
    l2(k_on) = k(k_on) - l1(k_on)

    # Model for cumulative localizations over time
    function model(x)
        n_emit, k_on_val = x
        n_emit_ceil = ceil(n_emit)
        k_val = k(k_on_val)
        l1_val = l1(k_on_val)
        l2_val = l2(k_on_val)

        return n_emit_ceil * (1.0 - p_miss) * (k_on_val / k_val) .*
               ((1 / l1_val) * (1.0 .- exp.(-l1_val .* (frames .- 1.0))) .-
                (1 / l2_val) * (1.0 .- exp.(-l2_val .* (frames .- 1.0))))
    end

    # Cost function (mean squared error)
    costfunction(x) = mean((nloccumulative .- model(x)).^2)

    # Optimization bounds
    lowerbound = [maximum(nlocperframe), 1e-5]
    upperbound = [Float64(nclusters), nloccumulative[end] / frames[end]]

    # Improved initial guess for number of emitters
    x1_guess = ceil(nclusters * k_bleach / (k_off * (1.0 - p_miss)))
    x1_guess = (x1_guess > lowerbound[1]) && (x1_guess < upperbound[1]) ?
               x1_guess : (upperbound[1] - lowerbound[1]) / 2.0

    initialguess = [x1_guess, 1.0 / frames[end]]

    # Optimize
    inner_optimizer = GradientDescent()
    optimresults = optimize(
        costfunction,
        lowerbound,
        upperbound,
        initialguess,
        Fminbox(inner_optimizer)
    )

    # Extract results
    xhat = optimresults.minimizer
    nemitters = xhat[1]
    k_on = xhat[2]

    return k_on, k_off, k_bleach, p_miss, nemitters
end
