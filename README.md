# SMLMFrameConnection

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMFrameConnection.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMFrameConnection.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMFrameConnection.jl)


## Overview
SMLMFrameConnection performs frame-connection on localization microscopy data organized in an SMLMData.SMLD2D structure (https://github.com/JuliaSMLM/SMLMData.jl). 
Specifically, SMLMFrameConnection connects repeated localizations of a single blinking event of an emitter into a single higher precision localization.  This is done 
using the algorithm(s) presented in https://doi.org/10.3389/fbinf.2021.724325.  Localizations which were connected will share the same unique integer value for the field
SMLMData.SMLD2D.connectID.


## Interface
Once an SMLMData.SMLD2D structure is fully populated, the user only needs to run a single high-level method from the package: `frameconnect()`.
All fields of the SMLMData.SMLD2D structure must be populated with either meaningful values (e.g., for fields like `x`, `y`, `σ_x`, `σ_y`, and `framenum`, 
which the algorithm depends on) or by placeholders with a meaningful size (e.g., fields like `bg` and `σ_bg`, which may not be available, should be set to something like 
`smld.bg = zeros(Float64, length(smld.framenum))`, and `σ_bg = fill(Inf64, length(smld.framenum))`).  

The algorithm can be run on the fully populated `smld::SMLMData.SMLD2D` with default parameters as

```smld_connected, smld_preclustered, smld_combined, params = SMLMFrameConnection.frameconnect(smld)```

The output `smld_connected` is a copy of `smld` with the field `smld.connectID` updated to associate connected localizations.  `smld_preclustered` is a copy of `smld` 
with the field `smld_preclustered.connectID` updated to associate localizations that belonged to the same precluster.  `smld_combined` contains the combined higher 
precision localizations (i.e., `smld_combined=SMLMFrameConnection.combinelocalizations(smld_connected)`) and is considered the main output of this package.  `params` 
is a structure containing the user-defined and internally-defined parameters used in the algorithm.

Several user-defined parameters can be changed as optional keyword arguments to `frameconnect()`:

```
smld_connected, smld_preclustered, smld_combined, params = SMLMFrameConnection.frameconnect(smld;
    nnearestclusters = 2, nsigmadev = 5.0,
    maxframegap = 5, nmaxnn = 2)
```
See the documentation in SMLMFrameConnection/structdefinitions.jl for a description of these parameters and for guidance to which functions use them.

## Citation
David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", 
Frontiers in Bioinformatics, 2021 
https://doi.org/10.3389/fbinf.2021.724325
