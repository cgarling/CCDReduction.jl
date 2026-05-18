```@meta
CurrentModule = CCDReduction
```

# Cosmic Ray Rejection

Cosmic rays produce sharp, extremely bright pixel hits in CCD images.
This example demonstrates how to detect and remove them using
[`cosmicray_lacosmic`](@ref).

## Synthetic image with cosmic rays

```@example cr
using CCDReduction
using Statistics

# Create a 100×100 sky background
background = fill(500.0, 100, 100)
# Add faint Gaussian noise
background .+= 5.0 .* randn(100, 100)

# Plant three cosmic ray hits of varying intensities
background[30, 40]  = 50_000.0
background[60, 70]  = 30_000.0
background[80, 20]  = 75_000.0

println("Max before cleaning: ", maximum(background))
println("Median background:   ", round(median(background), digits=1))
```

## Running LA Cosmic

```@example cr
cleaned, mask = cosmicray_lacosmic(background;
                                    sigma_clip = 4.5,
                                    niter      = 4,
                                    read_noise = 5.0,
                                    gain       = 1.0)

println("Cosmic rays detected: ", sum(mask))
println("Max after cleaning:   ", round(maximum(cleaned), digits=1))
println("Pixel (30,40) cleaned: ", round(cleaned[30, 40], digits=1))
println("Pixel (60,70) cleaned: ", round(cleaned[60, 70], digits=1))
println("Pixel (80,20) cleaned: ", round(cleaned[80, 20], digits=1))
```

## The cosmic ray mask

The second return value is a `BitMatrix`; `true` marks detected cosmic rays.

```@example cr
println("Flagged at (30,40): ", mask[30, 40])
println("Flagged at (60,70): ", mask[60, 70])
println("Flagged at (80,20): ", mask[80, 20])
# Show that an ordinary pixel is not flagged
println("Flagged at (10,10): ", mask[10, 10])
```

## Incorporating into a pipeline

[`cosmicray_lacosmic`](@ref) can be inserted directly into a
[`Pipeline`](@ref):

```@example cr
master_bias = fill(100.0, 100, 100)
master_flat = fill(  1.0, 100, 100)

pipe = Pipeline(
    f -> subtract_bias(f, master_bias),
    f -> flat_correct(f, master_flat),
    f -> first(cosmicray_lacosmic(f; gain=1.0, read_noise=5.0)),
)

result = process(pipe, background)
println("Result max: ", round(maximum(result), digits=1))
```

## Tuning the parameters

- **`sigma_clip`** – raise to be less aggressive (miss faint CRs), lower to
  flag more aggressively (risk false positives near bright compact sources).
- **`objlim`** – increase if real stars are being flagged; decrease if faint
  CRs are missed.
- **`niter`** – more iterations catch wider CR trails.
- **`read_noise`** and **`gain`** should always match your detector
  specifications for accurate noise estimation.
