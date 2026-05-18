# API Reference

Complete reference for all public functions and types in CCDReduction.jl.

## Index

```@index
```

## Data Types

```@docs
CCDData
CCDReduction.AbstractCCDData
CCDReduction.default_header
CCDReduction.CollectionRow
ImageCollection
CCDReduction.getdata
CCDReduction.writefits
```

## Calibration Methods

```@docs
subtract_bias
subtract_bias!
subtract_dark
subtract_dark!
subtract_overscan
subtract_overscan!
flat_correct
flat_correct!
gain_correct
gain_correct!
```

## Image Operations

```@docs
trim
trimview
crop
cropview
combine
```

## Advanced Reduction

```@docs
noise_model
cosmicray_lacosmic
```

## Pipelines

```@docs
Pipeline
ChunkedPipeline
tiled_step
process
```

## Collections

```@docs
fitscollection
arrays
filenames
ccds
Base.iterate(::ImageCollection)
Base.size(::ImageCollection)
Base.getindex(::ImageCollection, ::Integer)
Base.getindex(::ImageCollection, ::AbstractVector{Bool})
```

## Internal Utilities

```@docs
CCDReduction._laplacian
CCDReduction._median3x3
CCDReduction._median5x5
CCDReduction._dilate
CCDReduction.find_ccd
```
