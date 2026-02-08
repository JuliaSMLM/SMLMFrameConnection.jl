# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

SMLMFrameConnection performs frame-connection on 2D SMLM localization data: combining repeated localizations of a single blinking fluorophore across frames into a single higher-precision localization. Implements the spatiotemporal LAP algorithm from Schodt & Lidke (2021).

## Development Commands

```bash
# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file
julia --project=. -e 'using Test, SMLMFrameConnection, SMLMData; include("test/test_helpers.jl"); include("test/test_frameconnect.jl")'

# Build docs locally
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl

# Quick REPL usage
julia --project=.
```

Test files must `include("test/test_helpers.jl")` first since `runtests.jl` loads shared fixtures from it.

## Architecture

### Algorithm Pipeline

`frameconnect()` in `frameconnect.jl` is the main entry point. It orchestrates a 4-stage pipeline:

1. **Precluster** (`precluster.jl`): Spatiotemporal clustering using KDTree nearest-neighbor search. Groups localizations within `nsigmadev * mean(σ)` distance and `maxframegap` frames. Stores cluster assignments in emitter `track_id` fields.

2. **Estimate parameters** (`estimateparams.jl`, `estimatedensities.jl`): Estimates photophysics rates (k_on, k_off, k_bleach, p_miss) from precluster statistics. Uses Optim.jl NelderMead to fit cumulative localization counts. `estimatedensities.jl` estimates local emitter density per cluster using KDTree neighbor distances.

3. **Connect via LAP** (`connectlocalizations.jl` -> `create_costmatrix.jl` -> `solveLAP.jl` -> `linkclusters.jl`): For each multi-emitter precluster, builds a 2Nx2N cost matrix with connection/birth/death blocks using negative log-likelihoods from spatial separation, observation probability, and photophysics. Solves with Hungarian.jl. `linkclusters.jl` updates `connectID` from LAP assignments.

4. **Combine** (`combinelocalizations.jl`): MLE weighted mean using full 2x2 covariance (precision-weighted). Produces higher-precision output localizations.

### Data Flow

- Input/output: `BasicSMLD` from SMLMData.jl containing `Emitter2DFit` emitters
- Internal representation: `organizeclusters()` converts to `Vector{Matrix{Float32}}` where each matrix is a cluster with columns: `[x, y, σ_x, σ_y, σ_xy, frame, dataset, connectID, sortindex]`
- `ParamStruct`: mutable struct accumulating estimated parameters through the pipeline
- Return: tuple `(combined::BasicSMLD, info::FrameConnectInfo)`

### Key Types

- `FrameConnectConfig <: AbstractSMLMConfig`: User-facing config (keyword-constructible)
- `FrameConnectInfo{T} <: AbstractSMLMInfo`: Algorithm output metadata (track assignments, rates, timing)
- `ParamStruct`: Internal mutable state for pipeline parameters (not exported for direct use)

### Dual API Pattern

`frameconnect()` accepts both kwargs and a `FrameConnectConfig` struct. The kwargs form constructs the config internally. Both `FrameConnectConfig` and `FrameConnectInfo` inherit from SMLMData abstract types.

## Dependencies

- **SMLMData.jl** (v0.7): Provides `BasicSMLD`, `Emitter2DFit`, `IdealCamera`, abstract config/info types
- **Hungarian.jl**: LAP solver in `solveLAP.jl`
- **NearestNeighbors.jl**: KDTree for preclustering and density estimation
- **Optim.jl**: Parameter estimation (Fminbox + NelderMead)

## Conventions

- Positions and uncertainties in microns
- Frame numbers are 1-based integers
- `track_id=0` means unconnected; populated values are compressed to `1:n_tracks`
- Mutating functions use `!` suffix with non-mutating wrappers (e.g., `connectlocalizations!`/`connectlocalizations`)
- Emitter precision type `ET` is derived from emitter field values, not SMLD type parameter
