#---------------------------------------------------------------------------------------
# Pipeline
#---------------------------------------------------------------------------------------

"""
    Pipeline(steps...)

A sequential image-processing pipeline.

Each step must be a callable with signature `(frame) -> processed_frame`.
Steps are applied in order: the output of step *i* becomes the input of step
*i+1*.

Call [`process`](@ref) (or the pipeline itself as a function) to apply all
steps to a frame.

# Examples
```julia
master_bias = zeros(512, 512)
master_flat = ones(512, 512)

pipeline = Pipeline(
    frame -> subtract_bias(frame, master_bias),
    frame -> flat_correct(frame, master_flat),
)

science = ones(512, 512)
result  = process(pipeline, science)
```
"""
struct Pipeline
    steps::Vector{Any}
    Pipeline(steps...) = new(collect(Any, steps))
end

"""
    process(pipeline, frame)

Apply every step in `pipeline` to `frame` in order, returning the final result.

This is equivalent to calling the pipeline as a function: `pipeline(frame)`.
"""
function process(p::Pipeline, frame)
    result = frame
    for step in p.steps
        result = step(result)
    end
    return result
end

(p::Pipeline)(frame) = process(p, frame)

Base.length(p::Pipeline) = length(p.steps)

function Base.show(io::IO, p::Pipeline)
    print(io, "Pipeline with $(length(p.steps)) step(s)")
end


#---------------------------------------------------------------------------------------
# ChunkedPipeline
#---------------------------------------------------------------------------------------

"""
    ChunkedPipeline(steps...; chunksize = (256, 256))

A tile-based image-processing pipeline that improves cache efficiency by
applying *all* steps to a small chunk before moving on to the next chunk.

Unlike [`Pipeline`](@ref), which processes the entire image at each step,
`ChunkedPipeline` decomposes the image into tiles of size `chunksize` and runs
the full step sequence on each tile independently.  This avoids loading the
full image through each step and can dramatically improve L2/L3 cache hit rates
for large frames.

Each step must be a callable with signature `(chunk) -> processed_chunk`, where
`chunk` is a (possibly smaller) submatrix.  For calibration-frame-based steps
(e.g. bias, flat), the calibration frames must be sliced inside the closure:

```julia
master_bias = zeros(4096, 4096)
master_flat = ones(4096, 4096)

pipeline = ChunkedPipeline(
    chunk -> subtract_bias(chunk, master_bias[parentindices(chunk)...]),
    chunk -> flat_correct(chunk, master_flat[parentindices(chunk)...]);
    chunksize = (512, 512),
)
science = rand(4096, 4096)
result  = process(pipeline, science)
```

Alternatively, [`tiled_step`](@ref) automates the index slicing for you:

```julia
pipeline = ChunkedPipeline(
    tiled_step(subtract_bias, master_bias),
    tiled_step(flat_correct, master_flat);
    chunksize = (512, 512),
)
```

# Notes
- The output element type is always `float(eltype(frame))`.
- The order of tiles follows column-major (Julia) memory layout.

# See Also
[`tiled_step`](@ref), [`Pipeline`](@ref), [`process`](@ref)
"""
struct ChunkedPipeline
    steps     :: Vector{Any}
    chunksize :: Tuple{Int,Int}
    function ChunkedPipeline(steps...; chunksize::Tuple{Int,Int} = (256, 256))
        all(s -> s > 0, chunksize) || throw(ArgumentError("chunksize must be positive"))
        new(collect(Any, steps), chunksize)
    end
end

Base.length(p::ChunkedPipeline) = length(p.steps)

function Base.show(io::IO, p::ChunkedPipeline)
    print(io, "ChunkedPipeline with $(length(p.steps)) step(s), chunksize=$(p.chunksize)")
end

"""
    process(pipeline::ChunkedPipeline, frame)

Apply `pipeline` to `frame` tile by tile, returning the result as a new array.

For [`tiled_step`](@ref) steps, the calibration frames are automatically sliced
to the current tile region.  Non-`tiled_step` callables receive the output of
the previous step.
"""
function process(p::ChunkedPipeline, frame::AbstractMatrix)
    nrows, ncols = size(frame)
    cr, cc = p.chunksize
    F = float(eltype(frame))
    result = Array{F}(undef, nrows, ncols)

    for cj in 1:cc:ncols
        for ci in 1:cr:nrows
            ri = ci:min(ci + cr - 1, nrows)
            rj = cj:min(cj + cc - 1, ncols)
            processed = @view frame[ri, rj]
            for step in p.steps
                if step isa TiledStep
                    # Always use original (ri, rj) for calibration slicing
                    processed = step(processed, ri, rj)
                else
                    processed = step(processed)
                end
            end
            result[ri, rj] .= processed
        end
    end
    return result
end

(p::ChunkedPipeline)(frame) = process(p, frame)


#---------------------------------------------------------------------------------------
# tiled_step helper
#---------------------------------------------------------------------------------------

struct TiledStep{F, C <: Tuple}
    f    :: F
    cals :: C
end

"""
    tiled_step(f, calibrations...)

Wrap a calibration-based reduction function so that it automatically slices its
calibration frames to match the current tile in a [`ChunkedPipeline`](@ref).

`f` must accept a frame as its first argument followed by one calibration frame
per element of `calibrations`.  When used inside a `ChunkedPipeline`, the
pipeline automatically supplies the tile region indices so that each calibration
frame is sliced correctly.  The step can also be called directly on a
`SubArray` (outside a pipeline), in which case `parentindices` is used.

# Example
```julia
master_bias = zeros(4096, 4096)
master_flat = ones(4096, 4096)

pipeline = ChunkedPipeline(
    tiled_step(subtract_bias, master_bias),
    tiled_step(flat_correct, master_flat);
    chunksize = (512, 512),
)
result = process(pipeline, science_frame)
```
"""
tiled_step(f, cals...) = TiledStep(f, cals)

# Called by ChunkedPipeline with explicit (ri, rj) tile indices
function (s::TiledStep)(chunk::AbstractArray, ri, rj)
    sliced = map(c -> @view(c[ri, rj]), s.cals)
    return s.f(chunk, sliced...)
end

# Standalone call on a SubArray (e.g. in user code outside a pipeline)
function (s::TiledStep)(chunk::SubArray)
    idxs = parentindices(chunk)
    sliced = map(c -> @view(c[idxs...]), s.cals)
    return s.f(chunk, sliced...)
end

# Fallback for non-view chunks (e.g. plain Array)
function (s::TiledStep)(chunk::AbstractArray)
    return s.f(chunk, s.cals...)
end
