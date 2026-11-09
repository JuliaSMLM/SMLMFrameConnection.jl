# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

SMLMFrameConnection performs frame-connection on 2D SMLM localization data: combining repeated localizations of a single blinking fluorophore across frames into a single higher-precision localization. Implements the spatiotemporal LAP algorithm from Schodt & Lidke (2021).

## Development Commands

```bash
# Run tests (lab layout: Core = test/*.jl; Long = test/long/, slow, lab machine)
julia --project=. -e 'using Pkg; Pkg.test()'
GROUP=Everything julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file (each file loads its own packages and test_helpers.jl)
julia --project=test test/test_frameconnect.jl

# Build docs locally
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl

# Quick REPL usage
julia --project=.
```

`test/runtests.jl` is the lab's template runner (do not edit): each `test/*.jl` runs in its own module, so every test file starts with `using Test, SMLMFrameConnection, SMLMData` and `include("test_helpers.jl")`. Groups are declared in `test/test_groups.toml`.

## Architecture

### Algorithm Pipeline

`frameconnect()` in `frameconnect.jl` is the main entry point. It orchestrates a 5-stage pipeline:

`frameconnect` first sorts a copy of the input by `(dataset, frame, x, y)`, ties broken by every other field (`_order_key`), and runs every stage on it (calibration and combination included), so results do not depend on input order; only `info.connected` is permuted back to the caller's order. Input carrying `z` is refused (`_check_2d`). A precluster above 2000 localizations triggers a `@warn` (the LAP is O(N³) in precluster size).

1. **Precluster** (`precluster.jl`): Spatiotemporal clustering using KDTree nearest-neighbor search. Links each localization to its `max_neighbors` nearest neighbors within `max_sigma_dist * (mean σ_i + mean σ_j)` among the localizations of its frame and the `max_frame_gap` frames before; preclusters are the connected components of those links (union-find). Stores cluster assignments in emitter `track_id` fields.

2. **Estimate parameters** (`estimateparams.jl`, `estimatedensities.jl`): Estimates photophysics rates (k_on, k_off, k_bleach, p_miss) from precluster statistics. Uses Optim.jl NelderMead to fit cumulative localization counts. `estimatedensities.jl` estimates local emitter density per cluster using KDTree neighbor distances; a neighbor center closer than 1e-6 × the cluster's localization error counts as the same center and is skipped (a rounding-level separation would otherwise give a density of ~1e29 per µm² and change the linking).

3. **Connect via LAP** (`connectlocalizations.jl` -> `create_costmatrix.jl` -> `solveLAP.jl` -> `linkclusters.jl`): For each multi-emitter precluster, builds a 2Nx2N cost matrix with connection/birth/death blocks using negative log-likelihoods from spatial separation, observation probability, and photophysics. Solves with Hungarian.jl. `linkclusters.jl` updates `connectID` from LAP assignments.

4. **Calibrate** (optional, `calibration.jl`): When `config.calibration` is set, analyzes frame-to-frame jitter within connected tracks to estimate motion variance (`σ_motion²`) and CRLB scale factor (`k²`). The binned WLS fit is per axis `Var(Δx) = A + B·(σ_x1² + σ_x2²)`, so `A = 2σ_motion²` (stored as `CalibrationResult.A`) and `B = k²`. Pairs with chi² under the current fit (full pair covariance, σ_xy included; a trial A < 0 enters the cut as 0) above the chi²(2) 99.7% point are dropped, truncated bin means are divided by the truncation factor, bins are weighted by the fit-predicted variance, and the fit iterates to convergence (unbiased under the model; a fixed cut on the uncalibrated chi² biased k toward 1). Applies corrected uncertainties `Σ_corrected = σ_motion² I + k² Σ_CRLB` before combination. Bins (`_calibration_bins`) are equal-width, assigned by comparison with the edges; reported precisions whose range is at most `_CAL_EQUAL_REL` (1e-12) of the largest are one bin. The free fit is solved centered at the weighted mean of the bin centers (the uncentered normal equations cancel when the reported precisions barely vary). When it is undetermined or impossible (fewer than 3 bins), both one-parameter fits run (scale alone, motion alone); the scale alone is applied unless the motion term lowers the bins' chi-square (`_fit_chi2`: each bin mean's model variance from `_model_var_of_mean`), both scored on the pairs both cuts keep, by more than `_CAL_MOTION_CHI2` (3.84, the χ²(1) 95% point); with one bin both fit exactly and the scale is applied. Every applied fit must use at least 100 pairs (`_CAL_MIN_PAIRS`, counted in the bins the fit used, after the cut); otherwise it falls back (identity correction). A one-bin fit's standard error is model-based: each pair's covariance shape (`_pair_obs_var`) plus the correlation of consecutive pairs sharing a localization (`_model_var_of_mean`, computed in the fit loop only when a single bin survives, and once per candidate for `_fit_chi2`). Pair fields are promoted to Float64 before any arithmetic.

Covariance numerics (`structdefinitions.jl`): any positive-definite input stays positive definite through every step, however close |σ_xy| is to σ_xσ_y. `_cov_det` (Float64, exact sign via fma, factored) is the entry check and every single determinant; `_det_add`/`_pair_cov_det` give det of a sum as a sum of non-negative terms (pair costs in `create_costmatrix`, the calibration cut); `_precision_mean` (combine, density centers) carries the summed precision's determinant, takes x relative to the localization with the largest x precision and y relative to the one with the largest y precision, and gives identical covariances their exact closed form (x₁ + Σ(xᵢ - x₁)/n, Σ/n; a group of one comes out unchanged); `_round_cov` rounds a covariance to the emitter type keeping it positive definite (bounded steps). Chi² and Mahalanobis distances use the sum-of-squares form `dx²/sxx + sxx(dy - (sxy/sxx)dx)²/det`.

5. **Combine** (`combinelocalizations.jl`): MLE weighted mean using full 2x2 covariance (precision-weighted). Produces higher-precision output localizations. The `track_length` config (inclusive `(min, max)` range, `nothing`=disabled) drops tracks whose localization count falls outside the range here (offset-safe: `ncumulative` spans all localizations, a `keep` mask selects which clusters to emit) — `info.connected` retains all tracks, so filtering only affects the combined output. `info.n_filtered` reports the count dropped.

### Data Flow

- Input/output: `BasicSMLD` from SMLMData.jl containing `Emitter2DFit` emitters
- Internal representation: `organizeclusters()` converts to `Vector{Matrix{Float64}}` (Float64 so the integer columns stay exact above 2^24) where each matrix is a cluster, rows sorted by frame, with columns: `[x, y, σ_x, σ_y, σ_xy, frame, dataset, connectID, sortindex]`
- `ParamStruct`: mutable struct accumulating estimated parameters through the pipeline
- Return: tuple `(combined::BasicSMLD, info::FrameConnectInfo)`

### Key Types

- `FrameConnectConfig <: AbstractSMLMConfig`: User-facing config (keyword-constructible)
- `FrameConnectInfo{T} <: AbstractSMLMInfo`: Algorithm output metadata (track assignments, rates, timing)
- `CalibrationConfig`: Optional uncertainty calibration settings (nested in `FrameConnectConfig.calibration`)
- `CalibrationResult`: Calibration diagnostics (nested in `FrameConnectInfo.calibration`, `nothing` when disabled)
- `ParamStruct`: Internal mutable state for pipeline parameters (not exported; access via `SMLMFrameConnection.ParamStruct`)

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
