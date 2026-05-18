```@meta
CurrentModule = CCDReduction
```

# Basic Reduction Pipeline

This example walks through a complete CCD reduction sequence — from raw
overscan-containing frames to flat-field-corrected science images.  All data
are generated synthetically so the example can be run without any external
files.

## Synthetic data

```@example basic
using CCDReduction
using Statistics

# Image size: 100 rows × 110 columns (the last 10 columns are the overscan)
NROWS, NCOLS, NOVERSCAN = 100, 100, 10
SZ = (NROWS, NCOLS + NOVERSCAN)

# Read noise (σ) and gain (e⁻/ADU)
READ_NOISE = 10.0
GAIN       = 2.0

function make_raw(signal, bias_level = 500.0)
    shot_noise  = sqrt.(max.(signal ./ GAIN, 0.0)) .* randn(size(signal)...)
    read_noise  = READ_NOISE .* randn(size(signal)...)
    overscan    = fill(bias_level, NROWS, NOVERSCAN)
    # Combine image region + overscan
    hcat(signal .+ shot_noise .+ read_noise .+ bias_level, overscan)
end

# Master calibration frames (already reduced — normally these come from
# dedicated calibration exposures combined with `combine`)
master_bias = make_raw(zeros(NROWS, NCOLS)) |> x -> x[:, 1:NCOLS]
master_dark = make_raw(fill(5.0, NROWS, NCOLS); bias_level = 500.0) |> x -> x[:, 1:NCOLS]
master_flat = make_raw(fill(1000.0, NROWS, NCOLS); bias_level = 500.0) |> x -> x[:, 1:NCOLS]

# Three synthetic science frames with a faint source at (50,50) of signal 100
science_signal = fill(200.0, NROWS, NCOLS)
science_signal[48:52, 48:52] .= 500.0
raw_frames = [make_raw(science_signal) for _ in 1:3]

println("Raw frame size:  ", size(raw_frames[1]))
println("Calibration size: ", size(master_bias))
```

## Step 1 – Overscan subtraction & trimming

```@example basic
trimmed = [
    begin
        os_sub = subtract_overscan(f, (:, NCOLS+1:NCOLS+NOVERSCAN))
        trim(os_sub, (:, NCOLS+1:NCOLS+NOVERSCAN))
    end
    for f in raw_frames
]
println("Trimmed frame size: ", size(trimmed[1]))
```

## Step 2 – Bias subtraction

```@example basic
debiased = [subtract_bias(f, master_bias) for f in trimmed]
println("Mean after bias subtraction: ", round(mean(debiased[1]), digits=1))
```

## Step 3 – Dark subtraction

```@example basic
# Science exposure: 60 s; dark calibration exposure: 300 s
dark_corrected = [
    subtract_dark(f, master_dark; data_exposure=60.0, dark_exposure=300.0)
    for f in debiased
]
```

## Step 4 – Flat-field correction

```@example basic
flat_corrected = [flat_correct(f, master_flat) for f in dark_corrected]
println("Min: ", round(minimum(flat_corrected[1]), digits=2),
        "  Max: ", round(maximum(flat_corrected[1]), digits=2))
```

## Step 5 – Gain correction

```@example basic
electrons = [gain_correct(f, GAIN) for f in flat_corrected]
```

## Step 6 – Noise estimation

```@example basic
noise_maps = [noise_model(f; read_noise=READ_NOISE, gain=GAIN) for f in flat_corrected]
println("Typical noise (σ): ", round(mean(noise_maps[1]), digits=2), " ADU")
```

## Step 7 – Combine frames

```@example basic
combined = combine(flat_corrected)
println("Combined frame size: ", size(combined))
println("Combined (median) value near source: ",
        round(mean(combined[48:52, 48:52]), digits=1))
```

## Using Pipeline

The above can be expressed concisely as a [`Pipeline`](@ref):

```@example basic
pipe = Pipeline(
    f -> subtract_overscan(f, (:, NCOLS+1:NCOLS+NOVERSCAN)),
    f -> trim(f, (:, NCOLS+1:NCOLS+NOVERSCAN)),
    f -> subtract_bias(f, master_bias),
    f -> subtract_dark(f, master_dark; data_exposure=60.0, dark_exposure=300.0),
    f -> flat_correct(f, master_flat),
    f -> gain_correct(f, GAIN),
)

results = [process(pipe, f) for f in raw_frames]
combined2 = combine(results)
println("Pipeline result near source: ",
        round(mean(combined2[48:52, 48:52]), digits=1))
```
