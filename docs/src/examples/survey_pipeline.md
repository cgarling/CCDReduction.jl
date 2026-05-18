```@meta
CurrentModule = CCDReduction
```

# Batch Processing with ImageCollection

For a real survey you may have thousands of frames in a directory tree.
[`fitscollection`](@ref) scans all FITS files, reads their headers into an
[`ImageCollection`](@ref), and lets you filter and iterate over them with
minimal boilerplate.

This example generates a small mock survey on-the-fly.

## Creating mock survey data

```@example survey
using CCDReduction
using FITSIO
using Statistics

# Helper: write a 2-D array to a new single-HDU FITS file
function write_fits(path, data; kw...)
    hdr = FITSHeader(
        ["SIMPLE", string.(keys(kw))...],
        [true, values(kw)...],
        ["",   fill("", length(kw))...],
    )
    FITS(path, "w") do f
        write(f, permutedims(data, 2:-1:1); header = hdr)
    end
end

survey_dir = mktempdir()

# Generate 3 bias, 2 dark, 3 flat and 6 science frames
for i in 1:3
    write_fits(joinpath(survey_dir, "bias_$(i).fits"),
               fill(Float32(500), 64, 64);
               IMAGETYP = "BIAS", EXPTIME = 0.0)
end
for i in 1:2
    write_fits(joinpath(survey_dir, "dark_$(i).fits"),
               fill(Float32(2), 64, 64);
               IMAGETYP = "DARK", EXPTIME = 300.0)
end
for i in 1:3, filt in ("V", "R")
    write_fits(joinpath(survey_dir, "flat_$(filt)_$(i).fits"),
               fill(Float32(1000), 64, 64);
               IMAGETYP = "FLAT", FILTER = filt, EXPTIME = 30.0)
end
for i in 1:3, filt in ("V", "R")
    signal = fill(Float32(3000), 64, 64)
    signal[30:34, 30:34] .= Float32(8000)   # fake star
    write_fits(joinpath(survey_dir, "sci_$(filt)_$(i).fits"),
               signal;
               IMAGETYP = "LIGHT", FILTER = filt, EXPTIME = 120.0)
end

println("Files written: ", length(readdir(survey_dir)))
```

## Scanning the directory

```@example survey
col = fitscollection(survey_dir)
println(col)
```

## Filtering

```@example survey
bias_col    = col[col.IMAGETYP .== "BIAS"]
dark_col    = col[col.IMAGETYP .== "DARK"]
flat_V_col  = col[(col.IMAGETYP .== "FLAT")  .& (col.FILTER .== "V")]
flat_R_col  = col[(col.IMAGETYP .== "FLAT")  .& (col.FILTER .== "R")]
sci_V_col   = col[(col.IMAGETYP .== "LIGHT") .& (col.FILTER .== "V")]
sci_R_col   = col[(col.IMAGETYP .== "LIGHT") .& (col.FILTER .== "R")]

println("Bias frames:    ", length(bias_col))
println("Dark frames:    ", length(dark_col))
println("V-band flats:   ", length(flat_V_col))
println("R-band flats:   ", length(flat_R_col))
println("V-band science: ", length(sci_V_col))
println("R-band science: ", length(sci_R_col))
```

## Building master calibration frames

```@example survey
master_bias = combine(collect(arrays(bias_col)))
master_dark = combine(collect(arrays(dark_col)))
master_flat_V = combine(collect(arrays(flat_V_col)))
master_flat_R = combine(collect(arrays(flat_R_col)))

println("Master bias median: ", round(median(master_bias), digits=1))
```

## Processing science frames

Build a per-filter pipeline and apply it to all science frames at once:

```@example survey
function make_pipeline(flat, bias = master_bias, dark = master_dark)
    Pipeline(
        f -> subtract_bias(f, bias),
        f -> subtract_dark(f, dark; data_exposure=120.0, dark_exposure=300.0),
        f -> flat_correct(f, flat),
    )
end

pipe_V = make_pipeline(master_flat_V)
pipe_R = make_pipeline(master_flat_R)

# Process and collect results
reduced_V = [process(pipe_V, f) for f in arrays(sci_V_col)]
reduced_R = [process(pipe_R, f) for f in arrays(sci_R_col)]

combined_V = combine(reduced_V)
combined_R = combine(reduced_R)

println("V combined — star region mean: ",
        round(mean(combined_V[30:34, 30:34]), digits=2))
println("R combined — star region mean: ",
        round(mean(combined_R[30:34, 30:34]), digits=2))
```

## Saving results

[`ccds`](@ref) supports saving processed frames with a prefix/suffix:

```@example survey
out_dir = mktempdir()

ccds(sci_V_col; path = out_dir, save_prefix = "reduced_V") do img
    process(pipe_V, img)
end

out_col = fitscollection(out_dir; recursive = false)
println("Saved files: ", length(out_col))
println("First filename: ", out_col[1].name)
```
