```@meta
CurrentModule = CCDReduction
```

# CCDReduction.jl

CCDReduction.jl provides a complete toolkit for reducing optical imaging data
from CCD detectors, from single-frame calibration operations through full
survey-scale batch pipelines.

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/CCDReduction.jl)
[![Build Status](https://github.com/juliaastro/CCDReduction.jl/workflows/CI/badge.svg?branch=main)](https://github.com/juliaastro/CCDReduction.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/CCDReduction.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Codecov](https://codecov.io/gh/juliaastro/CCDReduction.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliaastro/CCDReduction.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation

From Julia enter Pkg mode (by pressing `]` in the Julia REPL)

```julia-repl
julia> ]

pkg> add CCDReduction
```

## Quick start

```@example quickstart
using CCDReduction

# Synthetic science + calibration frames
science = rand(Float64, 512, 512) .+ 500.0
master_bias = fill(100.0, 512, 512)
master_dark = fill(  2.0, 512, 512)   # 2 e⁻/s at 1-second dark exposure
master_flat = fill(  1.0, 512, 512)

# One-line reductions
debiased = subtract_bias(science, master_bias)
dark_sub = subtract_dark(debiased, master_dark; data_exposure=60.0, dark_exposure=1.0)
reduced  = flat_correct(dark_sub, master_flat)
nothing # hide
```

### Pipeline

Use [`Pipeline`](@ref) to chain steps and apply them to every frame:

```@example quickstart
pipe = Pipeline(
    f -> subtract_bias(f, master_bias),
    f -> subtract_dark(f, master_dark; data_exposure=60.0, dark_exposure=1.0),
    f -> flat_correct(f, master_flat),
)

result = process(pipe, science)
nothing # hide
```

For large frames, [`ChunkedPipeline`](@ref) applies all steps to one tile at a
time, improving cache efficiency:

```@example quickstart
chunked = ChunkedPipeline(
    tiled_step(subtract_bias, master_bias),
    tiled_step(subtract_dark, master_dark),
    tiled_step(flat_correct,  master_flat);
    chunksize = (128, 128),
)

result2 = process(chunked, science)
nothing # hide
```

### Batch processing with ImageCollection

Scan a directory of FITS files and iterate over them easily:

```julia
col = fitscollection("data/")

# access columns
col.paths
col.EXPTIME         # all EXPTIME values as a Vector

# filter to science frames
sci = col[col.IMAGETYP .== "LIGHT"]

# process and save
ccds(sci; path = "reduced/", save_prefix = "r") do img
    pipe(img)
end
```

## Features

| Feature | Function |
|---------|----------|
| Bias subtraction | [`subtract_bias`](@ref) |
| Overscan subtraction | [`subtract_overscan`](@ref) |
| Dark subtraction | [`subtract_dark`](@ref) |
| Flat-field correction | [`flat_correct`](@ref) |
| Gain correction | [`gain_correct`](@ref) |
| Frame trimming | [`trim`](@ref) / [`trimview`](@ref) |
| Frame cropping | [`crop`](@ref) / [`cropview`](@ref) |
| Frame combination | [`combine`](@ref) |
| Noise estimation | [`noise_model`](@ref) |
| Cosmic-ray rejection | [`cosmicray_lacosmic`](@ref) |
| Batch collection | [`fitscollection`](@ref) / [`ImageCollection`](@ref) |
| Sequential pipeline | [`Pipeline`](@ref) |
| Cache-friendly pipeline | [`ChunkedPipeline`](@ref) |

## License

This work is distributed under the MIT "expat" license. See
[`LICENSE`](https://github.com/juliaastro/CCDReduction.jl/blob/main/LICENSE)
for more information.
