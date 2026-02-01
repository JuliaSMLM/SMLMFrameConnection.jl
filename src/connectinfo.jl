"""
    ConnectInfo{T}

Secondary output from `frameconnect()` containing connected localizations and metadata.

# Fields
- `connected::BasicSMLD{T}`: Input SMLD with `track_id` assigned (localizations uncombined)
- `n_input::Int`: Number of input localizations
- `n_tracks::Int`: Number of tracks formed
- `n_combined::Int`: Number of output (combined) localizations
- `k_on::Float64`: Estimated on rate (1/frame)
- `k_off::Float64`: Estimated off rate (1/frame)
- `k_bleach::Float64`: Estimated bleach rate (1/frame)
- `p_miss::Float64`: Probability of missed detection
- `initialdensity::Vector{Float64}`: Initial density estimate per precluster (emitters/μm²)
- `elapsed_ns::UInt64`: Wall time in nanoseconds
- `algorithm::Symbol`: Algorithm used (`:lap`)
- `n_preclusters::Int`: Number of preclusters formed
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
    elapsed_ns::UInt64
    algorithm::Symbol
    n_preclusters::Int
end
