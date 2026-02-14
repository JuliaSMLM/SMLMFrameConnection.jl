"""
    costmatrix = create_costmatrix(clusterdata::Vector{Matrix{Float32}},
        params::ParamStruct, clusterind::Int64, nframes::Int64,
        base_rho_on::Vector{Float64})

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

`base_rho_on` is a precomputed unit-density on-state density curve of length
`nframes`, shared across all clusters to avoid per-cluster allocation.
"""
function create_costmatrix(clusterdata::Vector{Matrix{Float32}},
    params::ParamStruct, clusterind::Int64, nframes::Int64,
    base_rho_on::Vector{Float64})
    # Extract some entries from clusterdata and params to make the code more
    # readable.
    # Column layout: 1:x 2:y 3:σ_x 4:σ_y 5:σ_xy 6:frame 7:dataset 8:connectID 9:sortindex
    cluster = clusterdata[clusterind]
    x = @view cluster[:, 1]
    y = @view cluster[:, 2]
    x_se = @view cluster[:, 3]
    y_se = @view cluster[:, 4]
    xy_cov = @view cluster[:, 5]
    framenum = @view cluster[:, 6]
    k_on = params.k_on
    k_off = params.k_off
    k_bleach = params.k_bleach
    p_miss = params.p_miss
    rho_0 = params.initial_density[clusterind]
    max_frame_gap = params.max_frame_gap

    # Populate the upper-left "connection" block.
    nlocalizations = length(framenum)
    n2 = 2 * nlocalizations
    costmatrix = fill(Inf32, (n2, n2))
    validcost = fill(false, (n2, n2))
    for mm = 1:nlocalizations, nn = (mm+1):nlocalizations
        deltaframe = abs(framenum[mm] - framenum[nn])
        if deltaframe == 0
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
    # Use precomputed base_rho_on, scaled by per-cluster density rho_0.
    startframe = 1
    birth_range = (nlocalizations+1):n2
    for nn in 1:nlocalizations
        frm = Int32(framenum[nn])
        deltaframe_past = min(max_frame_gap, frm - startframe)
        deltaframe_future = min(max_frame_gap, nframes - frm)
        # Localization uncertainty ellipse area (μm²) to make density dimensionless
        # Area = π√det(Σ) where det(Σ) = σ_x²σ_y² - σ_xy²
        det_loc = x_se[nn]^2 * y_se[nn]^2 - xy_cov[nn]^2
        loc_area = pi * sqrt(max(det_loc, eps()))
        rho_off_f = rho_0 * base_rho_on[frm] * k_off / k_on
        rho_on_past = rho_0 * base_rho_on[frm - deltaframe_past]
        birthcost = -log(1-p_miss) -
            log(rho_off_f * loc_area * (1-exp(-k_on)) * exp(-deltaframe_past*k_on) +
                rho_on_past * loc_area * (p_miss^deltaframe_past))
        deathcost = -log((1-exp(-k_off)) +
            (1-exp(-k_bleach)) +
            (p_miss^deltaframe_future))
        costmatrix[birth_range, nn] .= birthcost
        costmatrix[nn, birth_range] .= deathcost
        validcost[birth_range, nn] .= true
        validcost[nn, birth_range] .= true
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
