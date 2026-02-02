# ConnectInfo struct for tuple-pattern return

"""
    ConnectInfo{T}

Secondary output from `frameconnect()` containing track assignments and algorithm metadata.

# Fields
- `connected::BasicSMLD{T}`: Input SMLD with track_id assigned (localizations uncombined)
- `n_input::Int`: Number of input localizations
- `n_tracks::Int`: Number of tracks formed
- `n_combined::Int`: Number of output localizations
- `k_on::Float64`: Estimated on rate (1/frame)
- `k_off::Float64`: Estimated off rate (1/frame)
- `k_bleach::Float64`: Estimated bleach rate (1/frame)
- `p_miss::Float64`: Probability of missed detection
- `initialdensity::Vector{Float64}`: Initial density estimate per cluster (emitters/μm²)
- `elapsed_s::Float64`: Wall time in seconds
- `algorithm::Symbol`: Algorithm used (`:lap`)
- `n_preclusters::Int`: Number of preclusters formed

# Rate Parameter Interpretation

The rate parameters describe the photophysics of blinking fluorophores:
- `k_on`: Rate at which dark emitters convert to visible state
- `k_off`: Rate at which visible emitters convert to reversible dark state
- `k_bleach`: Rate at which visible emitters are irreversibly photobleached
- Duty cycle = k_on / (k_on + k_off)
- For typical dSTORM: k_on << k_off (low duty cycle, mostly dark with brief blinks)

# Example
```julia
(combined, info) = frameconnect(smld)
println("Connected \$(info.n_input) → \$(info.n_combined) localizations")
println("Formed \$(info.n_tracks) tracks in \$(info.elapsed_s)s")
# Access track assignments via info.connected
```
"""
struct ConnectInfo{T}
    connected::BasicSMLD{T}
    n_input::Int
    n_tracks::Int
    n_combined::Int
    k_on::Float64
    k_off::Float64
    k_bleach::Float64
    p_miss::Float64
    initialdensity::Vector{Float64}
    elapsed_s::Float64
    algorithm::Symbol
    n_preclusters::Int
end
