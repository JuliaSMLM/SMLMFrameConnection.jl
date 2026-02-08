"""
    costmatrix = create_costmatrix(clusterdata::Vector{Matrix{Float32}}, 
        params::ParamStruct, nframes::Int64)

Create a cost matrix for connections between localizations in `clusterdata`.

# Description
The costs associated with connecting/not connecting the localizations in the
input set of data in `clusterdata` are placed in entries of the output matrix
`costmatrix`.  For N localizations in `clusterdata`, the output `costmatrix`
will be an 2Nx2N block matrix.  The NxN blocks are attributed the following
meanings: the upper-left "connection" block corresponds to adding localizations
to existing clusters of other localizations, the bottom-left "birth" block 
corresponds to introduction of a new cluster, and the upper-right "death" block
corresponds to preventing admission of any more localizations to a cluster.
The bottom-right "auxillary" block is the transpose of the "connection" block.
Forbidden connections are indicated by `missing`.
"""
function create_costmatrix(clusterdata::Vector{Matrix{Float32}}, 
    params::ParamStruct, clusterind::Int64, nframes::Int64)
    # Extract some entries from clusterdata and params to make the code more
    # readable.
    # Column layout: 1:x 2:y 3:σ_x 4:σ_y 5:σ_xy 6:frame 7:dataset 8:connectID 9:sortindex
    x = clusterdata[clusterind][:, 1]
    y = clusterdata[clusterind][:, 2]
    x_se = clusterdata[clusterind][:, 3]
    y_se = clusterdata[clusterind][:, 4]
    xy_cov = clusterdata[clusterind][:, 5]
    framenum = clusterdata[clusterind][:, 6]
    k_on = params.k_on
    k_off = params.k_off
    k_bleach = params.k_bleach
    p_miss = params.p_miss
    rho_0 = params.initial_density[clusterind]
    max_frame_gap = params.max_frame_gap

    # Populate the upper-left "connection" block.
    nlocalizations = length(framenum)
    costmatrix = fill(Inf32, (2*nlocalizations, 2*nlocalizations))
    validcost = fill(false, (2*nlocalizations, 2*nlocalizations))
    for mm = 1:nlocalizations, nn = (mm+1):nlocalizations
        # For each localization, define the cost of linking to the other
        # localizations (divided by two, since a selection of one of these
        # costs leads to the same selection in the auxillary block).  In this
        # case, it's the negative log-likelihood, where the likelihood is given
        # as p(observed separation | localization error) ...
        #  * p(missing last N frames localizations)p(not missing new loc.)...
        #       *p(not turning off or bleaching)
        # Costs are divided by two and split among the "connection" and
        # "auxillary" blocks.
        deltaframe = abs(framenum[mm] - framenum[nn])
        if deltaframe == 0
            # Localizations in the same frame should never be connected.
            # NOTE: preclustering allows localizations in the same frame to be
            #       part of the same cluster, which is why this block is
            #       needed.
            continue
        end
        # Combined covariance matrix for two localizations:
        # Σ_comb = Σ_mm + Σ_nn (variances add for independent measurements)
        σ_x² = x_se[mm]^2 + x_se[nn]^2
        σ_y² = y_se[mm]^2 + y_se[nn]^2
        σ_xy = xy_cov[mm] + xy_cov[nn]
        det_Σ = σ_x² * σ_y² - σ_xy^2

        # Mahalanobis distance: Δᵀ Σ⁻¹ Δ where Σ⁻¹ = [σ_y² -σ_xy; -σ_xy σ_x²] / det
        Δx, Δy = x[mm] - x[nn], y[mm] - y[nn]
        mahal = (σ_y² * Δx^2 - 2*σ_xy*Δx*Δy + σ_x² * Δy^2) / det_Σ

        # Negative log-likelihood: -log(p) = 0.5*log(2π) + 0.5*log(det) + 0.5*mahal
        separationcost = log(2*pi) + 0.5*log(det_Σ) + 0.5*mahal
        observationcost = -log((p_miss^(deltaframe-1)) * (1-p_miss))
        stilloncost = (k_off+k_bleach) * deltaframe
        costmatrix[mm, nn] = (separationcost+observationcost+stilloncost) / 2.0
        costmatrix[nn+nlocalizations, mm+nlocalizations] = costmatrix[mm, nn]
        validcost[mm, nn] = true
        validcost[nn+nlocalizations, mm+nlocalizations] = true
    end

    # Populate the lower-left "birth" block and the upper-right "death" block.
    startframe = 1
    dutycycle = k_on / (k_on+k_off+k_bleach)
    frames = collect(1:maximum(framenum))
    lambda1 = k_bleach * dutycycle
    lambda2 = (k_on+k_off+k_bleach) - lambda1
    rho_on = rho_0 .* dutycycle .* 
        (exp.(-lambda1*(frames.-1)) - exp.(-lambda2*(frames.-1)));
    rho_off = rho_on .* k_off / k_on
    indices = 1:nlocalizations
    framesint = Int32.(framenum)
    for nn in indices
        deltaframe_past = minimum([max_frame_gap; framesint[nn]-startframe])
        deltaframe_future = minimum([max_frame_gap; nframes-framesint[nn]])
        # Localization uncertainty ellipse area (μm²) to make density dimensionless
        # Area = π√det(Σ) where det(Σ) = σ_x²σ_y² - σ_xy²
        det_loc = x_se[nn]^2 * y_se[nn]^2 - xy_cov[nn]^2
        loc_area = pi * sqrt(max(det_loc, eps()))
        birthcost = -log(1-p_miss) -
            log(rho_off[framesint[nn]] * loc_area * (1-exp(-k_on)) * exp(-deltaframe_past*k_on) +
                rho_on[framesint[nn]-deltaframe_past] * loc_area * (p_miss^deltaframe_past))
        deathcost = -log((1-exp(-k_off)) +
            (1-exp(-k_bleach)) +
            (p_miss^deltaframe_future))
        costmatrix[indices .+ nlocalizations, nn] .= birthcost
        costmatrix[nn, indices .+ nlocalizations] .= deathcost
        validcost[indices .+ nlocalizations, nn] .= true
        validcost[nn, indices .+ nlocalizations] .= true
    end

    # If any of the valid assignments (```validcost```) are inf, we should
    # replace them with a usable value.
    badcosts = isinf.(costmatrix)
    if any(badcosts[validcost[:]])
        costsum = sum(costmatrix[.!badcosts])
        costmatrix[badcosts[:] .& validcost[:]] .= 2 * costsum
    end

    return costmatrix
end