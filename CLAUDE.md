# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SMLMFrameConnection.jl performs frame-connection on 2D localization microscopy data, combining repeated localizations of blinking fluorophores into higher-precision localizations. Part of the JuliaSMLM ecosystem.

**Reference:** Schodt & Lidke (2021), "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem"

## Commands

```bash
# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Run a specific test file
julia --project test/test_frameconnect.jl

# Interactive development
julia --project
using SMLMFrameConnection, SMLMData
```

## Architecture

### Algorithm Pipeline

The main entry point `frameconnect()` executes four stages:

1. **Preclustering** (`precluster.jl`): Groups spatiotemporally nearby localizations using KD-tree nearest-neighbor search. Distance threshold: `nsigmadev × mean(σ_x, σ_y)`. Temporal threshold: `maxframegap` frames.

2. **Parameter Estimation** (`estimateparams.jl`, `estimatedensities.jl`): Estimates photophysics (k_on, k_off, k_bleach, p_miss) and emitter densities from precluster statistics.

3. **LAP Solving** (`create_costmatrix.jl`, `solveLAP.jl`, `linkclusters.jl`): Constructs cost matrices for each precluster and solves linear assignment problems via Hungarian algorithm to determine which localizations belong to the same physical emitter.

4. **Combination** (`combinelocalizations.jl`): Merges connected localizations using MLE weighted mean: `x = Σ(x/σ²) / Σ(1/σ²)`, uncertainty: `σ = √(1/Σ(1/σ²))`.

### Key Types

- **Input**: `BasicSMLD{T, Emitter2DFit{T}}` from SMLMData.jl
- **Output**: `(combined::BasicSMLD, info::ConnectInfo)` tuple
- **ConnectInfo{T}** (`connectinfo.jl`): Secondary output with connected SMLD, statistics, and photophysics
- **ParamStruct** (`structdefinitions.jl`): Internal algorithm parameters

### File Organization

```
src/
├── frameconnect.jl          # Main entry point
├── connectinfo.jl           # ConnectInfo output struct
├── structdefinitions.jl     # ParamStruct, helper conversions
├── precluster.jl            # Stage 1: spatiotemporal clustering
├── organizeclusters.jl      # Prepare cluster data structures
├── estimateparams.jl        # Stage 2: photophysics estimation
├── estimatedensities.jl     # Stage 2: density estimation
├── create_costmatrix.jl     # Stage 3: LAP cost matrix construction
├── solveLAP.jl              # Stage 3: Hungarian algorithm wrapper
├── linkclusters.jl          # Stage 3: apply LAP solution
├── combinelocalizations.jl  # Stage 4: MLE combination with covariance
├── defineidealFC.jl         # Ground-truth validation utility
└── compress_connectID.jl    # Normalize track_id to 1:N
```

### Dependencies

- **SMLMData.jl** (v0.6+): Core types (BasicSMLD, Emitter2DFit with σ_xy)
- **Hungarian.jl**: LAP solver
- **NearestNeighbors.jl**: KD-tree for spatial queries
- **Optim.jl**: Parameter optimization

## Testing

Tests use helper fixtures in `test/test_helpers.jl`:
- `make_test_camera()`, `make_emitter()`, `make_test_smld()`
- `make_blinking_molecule()`: Creates emitter sequence for one molecule
- `make_two_molecules_smld()`, `make_single_molecule_smld()`

## API

Five exports:
- `frameconnect(smld)` → `(combined::BasicSMLD, info::ConnectInfo)` tuple
- `ConnectInfo{T}` → struct with connected SMLD, statistics (n_input, n_tracks, n_combined), photophysics (k_on, k_off, k_bleach, p_miss), timing
- `combinelocalizations(smld)` → combines emitters by `track_id` using MLE with full covariance
- `defineidealFC(smld)` → ground-truth frame connection for simulated data
- `ParamStruct` → internal algorithm parameters
