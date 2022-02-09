using Hungarian

"""
    assignment, cost = solveLAP(costmatrix::Matrix{Float32})

Solve the linear assignment problem designated by `costmatrix`.

# Description
Solve the linear assignment problem (LAP) defined by `costmatrix` using the
package Hungarian.jl.  For now, this method is just a wrapper around 
Hungarian.jl, however the intention is to add more checks and options to this
method in the future.
"""
function solveLAP(costmatrix::Matrix{Float32})
    # If the input costmatrix is empty we shouldn't continue.
    if isempty(costmatrix)
        assignment = []
        cost = []
        return assignment, cost
    end

    # Solve the LAP.
    assignment, cost = Hungarian.hungarian(costmatrix)

    return assignment, cost
end