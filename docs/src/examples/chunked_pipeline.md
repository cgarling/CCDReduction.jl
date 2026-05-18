```@meta
CurrentModule = CCDReduction
```

# ChunkedPipeline and Cache Efficiency

For very large frames (e.g. 4k×4k or 8k×8k CCDs), loading the entire frame
into each pipeline step sequentially means every pixel is read and written
multiple times to RAM.  [`ChunkedPipeline`](@ref) improves cache efficiency by
processing one small tile at a time: *all* pipeline steps are applied to a
tile before the next tile is loaded.

## Concept

Consider three steps: bias subtraction → dark subtraction → flat correction.

**Sequential Pipeline** (loads the full frame at every step):
```
frame  →  step1(full)  →  step2(full)  →  step3(full)  →  result
```

**ChunkedPipeline** (processes one 256×256 tile at a time):
```
tile₁  →  step1  →  step2  →  step3  →  result[tile₁]
tile₂  →  step1  →  step2  →  step3  →  result[tile₂]
...
```

## Example

```@example chunked
using CCDReduction
using Statistics

# Synthetic 1024×1024 science + calibration frames
SZ = (1024, 1024)
science     = fill(2000.0, SZ...)
master_bias = fill( 500.0, SZ...)
master_dark = fill(   2.0, SZ...)
master_flat = fill(   1.0, SZ...)

# --- Sequential Pipeline ---------------------------------------------------
seq_pipe = Pipeline(
    f -> subtract_bias(f, master_bias),
    f -> subtract_dark(f, master_dark; data_exposure=120.0, dark_exposure=300.0),
    f -> flat_correct(f, master_flat),
)

# --- ChunkedPipeline (256×256 tiles) ----------------------------------------
chunk_pipe = ChunkedPipeline(
    tiled_step(subtract_bias, master_bias),
    tiled_step(subtract_dark, master_dark),
    tiled_step(flat_correct,  master_flat);
    chunksize = (256, 256),
)

seq_result   = process(seq_pipe,   science)
chunk_result = process(chunk_pipe, science)

println("Results identical: ", seq_result ≈ chunk_result)
println("Mean result: ", round(mean(chunk_result), digits=3))
```

## Mixing tiled_step with plain closures

[`tiled_step`](@ref) is only needed for steps that use a full-size calibration
frame.  Steps that apply a per-pixel operation without reference frames can be
written as plain closures and they interoperate seamlessly:

```@example chunked
gain     = 2.0
rn       = 8.0

mixed_pipe = ChunkedPipeline(
    tiled_step(subtract_bias, master_bias),
    tiled_step(flat_correct,  master_flat),
    chunk -> gain_correct(chunk, gain);   # plain closure, no calibration frame
    chunksize = (128, 128),
)

result2 = process(mixed_pipe, science)
println("Mean after gain correct: ", round(mean(result2), digits=1))
```

## Choosing chunksize

- Smaller chunks → better cache locality but more function-call overhead.
- A good starting point is 256×256 (262,144 `Float64` elements ≈ 2 MB),
  which fits comfortably in L2 cache on most modern CPUs.
- Experiment with the [`BenchmarkTools`](https://github.com/JuliaCI/BenchmarkTools.jl)
  package to find the optimal size for your hardware.
