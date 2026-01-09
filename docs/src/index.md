```@meta
CurrentModule = SMLMFrameConnection
```

# SMLMFrameConnection

Frame-connection for single molecule localization microscopy (SMLM) data. Combines repeated localizations of blinking events into higher-precision localizations using the algorithm from [Schodt & Lidke 2021](https://doi.org/10.3389/fbinf.2021.724325).

## Installation

```julia
using Pkg
Pkg.add("SMLMFrameConnection")
```

## Quick Start

```julia
using SMLMData, SMLMFrameConnection

# Create or load your BasicSMLD with Emitter2DFit emitters
smld_connected, smld_preclustered, smld_combined, params = frameconnect(smld)
```

The main output is `smld_combined`, which contains combined higher-precision localizations.

## Input Requirements

The input `BasicSMLD` must contain `Emitter2DFit` emitters with:
- `x`, `y`: Position coordinates (microns)
- `σ_x`, `σ_y`: Position uncertainties (microns) - required for MLE combination
- `frame`: Frame number
- `dataset`: Dataset identifier

## Algorithm Overview

1. **Preclustering**: Groups localizations based on spatial proximity and temporal adjacency
2. **Parameter Estimation**: Estimates blinking kinetics (k_on, k_off, k_bleach) and emitter density
3. **Frame Connection**: Uses Linear Assignment Problem (LAP) to optimally connect localizations
4. **Combination**: Combines connected localizations using MLE weighted mean

## API Reference

```@index
```

```@autodocs
Modules = [SMLMFrameConnection]
```

## Citation

David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Frontiers in Bioinformatics, 2021. [DOI: 10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)
