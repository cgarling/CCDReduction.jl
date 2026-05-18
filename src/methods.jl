# helper function
axes_min_length(idxs) = argmin([a isa Colon ? Inf : length(a) for a in idxs])

# parses FITS indices into standard Julian indices
function fits_indices(string::String)
    str = replace(string, r"[\[\]\s]" => "")
    tokens = split(str, ',')

    idxs = map(tokens) do token
        t = split(token, ':', keepempty=false)
        length(t) == 0 ? Colon() : parse(Int, t[1]):parse(Int, t[2])
    end
    return reverse(idxs)
end

# Convert using `round` for integers
convert_value(S::Type{<:Integer}, x) = round(S, x)
convert_value(S, x) = convert(S, x)

#-------------------------------------------------------------------------------
"""
    subtract_bias!(frame, bias_frame)

In-place version of [`subtract_bias`](@ref).

# See Also
[`subtract_bias`](@ref)
"""
function subtract_bias!(frame::AbstractArray, bias_frame::AbstractArray)
    frame .-= bias_frame
    return frame
end


"""
    subtract_bias(frame, bias_frame)

Subtract the `bias_frame` from `frame`.

If either are strings, they will be loaded into [`CCDData`](@ref) first. The HDU
loaded can be specified by `hdu` as either an integer or a tuple corresponding
to each file.

# Examples
```jldoctest
julia> frame = [1.0 2.2 3.3 4.5];

julia> bias = [0.0 0.2 0.3 0.5];

julia> subtract_bias(frame, bias)
1×4 Matrix{Float64}:
 1.0  2.0  3.0  4.0

```

# See Also
[`subtract_bias!`](@ref)
"""
subtract_bias(frame::AbstractArray, bias_frame::AbstractArray) = subtract_bias!(deepcopy(frame), bias_frame)


"""
    subtract_overscan!(frame, idxs; dims = axes_min_length(idxs))

In-place version of [`subtract_overscan`](@ref).

# See Also
[`subtract_overscan`](@ref)
"""
function subtract_overscan!(frame::AbstractArray{T}, idxs; dims = axes_min_length(idxs)) where T
    overscan_region = @view frame[idxs...]
    overscan_value = convert_value.(T, median(overscan_region, dims = dims))
    frame .-= overscan_value
    return frame
end

subtract_overscan!(frame::AbstractArray, idxs::String; kwargs...) = subtract_overscan!(frame, fits_indices(idxs); kwargs...)


"""
    subtract_overscan(frame, idxs; dims = axes_min_length(idxs))

Subtract the overscan frame from image.

`dims` is the dimension along which `overscan_frame` is combined. The default
value of `dims` is the axis with smaller length in overscan region.
If `idxs` is a string it will be parsed as FITS-style indices.

If `frame` is a string, it will be loaded into [`CCDData`](@ref) first.
The HDU loaded can be specified by `hdu` which by default is 1.

# Examples
```jldoctest
julia> frame = [4.0 2.0 3.0 1.0 1.0];

julia> subtract_overscan(frame, (:, 4:5), dims = 2)
1×5 Matrix{Float64}:
 3.0  1.0  2.0  0.0  0.0

julia> subtract_overscan(frame, "[4:5, 1:1]", dims = 2)
1×5 Matrix{Float64}:
 3.0  1.0  2.0  0.0  0.0
```

# See Also
[`subtract_overscan!`](@ref)
"""
subtract_overscan(frame, idxs; kwargs...) = subtract_overscan!(deepcopy(frame), idxs; kwargs...)


"""
    flat_correct!(frame, flat_frame; norm_value = mean(flat_frame))

In-place version of [`flat_correct`](@ref).

# See Also
[`flat_correct`](@ref)
"""
function flat_correct!(frame::AbstractArray, flat_frame::AbstractArray; norm_value = mean(flat_frame))
    norm_value <= 0 && error("norm_value must be positive")
    frame ./= (flat_frame ./ norm_value)
    return frame
end


"""
    flat_correct(frame, flat_frame; norm_value = mean(flat_frame))

Correct `frame` for non-uniformity using the calibrated `flat_frame`.

By default, the `flat_frame` is normalized by its mean, but this can be changed
by providing a custom `norm_value`.

If either are strings, they will be loaded into [`CCDData`](@ref) first. The HDU
loaded can be specified by `hdu` as either an integer or a tuple corresponding
to each file.

!!! note
    This function may introduce non-finite values if `flat_frame` contains
    values very close to `0` due to dividing by zero.
    The default behavior will return `Inf` if the frame value is non-zero, and
    `Nan` if the frame value is `0`.

# Examples
```jldoctest
julia> frame = ones(3, 3);

julia> flat = fill(2.0, (3, 3));

julia> flat_correct(frame, flat, norm_value = 1.0)
3×3 Matrix{Float64}:
 0.5  0.5  0.5
 0.5  0.5  0.5
 0.5  0.5  0.5

julia> flat_correct(frame, flat)
3×3 Matrix{Float64}:
 1.0  1.0  1.0
 1.0  1.0  1.0
 1.0  1.0  1.0
```

# See Also
[`flat_correct!`](@ref)
"""
function flat_correct(frame::AbstractArray{T}, flat_frame::AbstractArray{S}; kwargs...) where {T, S}
    V = float(promote_type(T, S))
    return flat_correct!(V.(frame), flat_frame; kwargs...)
end


"""
    trim(frame, idxs)

Trim the `frame` to remove the region specified by `idxs`.

This function trims the array in a manner such that the final array is
rectangular. The indices follow standard Julia convention, so `(:, 45:60)`
trims all columns from 45 to 60 and `(1:20, :)` trims all the rows from 1 to
20. The function also supports FITS-style index strings.

If `frame` is a string, it will be loaded into [`CCDData`](@ref) first.
The HDU loaded can be specified by `hdu` which by default is 1.

# Examples
```jldoctest
julia> frame = ones(5, 5);

julia> trim(frame, (:, 2:5))
5×1 Matrix{Float64}:
 1.0
 1.0
 1.0
 1.0
 1.0

julia> trim(frame, "[2:5, 1:5]")
5×1 Matrix{Float64}:
 1.0
 1.0
 1.0
 1.0
 1.0

```

# See Also
[`trimview`](@ref)
"""
trim(frame, idxs) = copy(trimview(frame, idxs))


"""
    trimview(frame, idxs)

Return a view of `frame` with the region specified by `idxs` removed.

This function is the same as [`trim`](@ref) but returns a view of the frame
rather than a copy.

!!! note
    Because this returns a view, any modification to the output will also
    modify `frame`.

# See Also
[`trim`](@ref)
"""
function trimview(frame::AbstractArray, idxs)
    # this adds the support for input indices of the form (1:size(frame, 1), ...) or (..., 1:size(frame, 2))
    # It converts 1:size(frame, 1) to : and then the same subroutine follows.
    processed_idxs = map(axes(frame), idxs) do a1, a2
        (a2 isa Colon || a1 == a2) ? Colon() : a2
    end

    # can switch to using `only` for Julia v1.4+
    ds = findall(x -> !isa(x, Colon), processed_idxs)
    length(ds) == 1 || error("Invalid trim indices $idxs")

    d = ds[1]
    full_idxs = axes(frame, d)
    # checking bounds error
    processed_idxs[d] ⊆ full_idxs || error("Trim indices $(idxs[d]) out of bounds for frame dimension $d $(full_idxs)")

    # finding the complement indices
    complement_idxs = setdiff(full_idxs, processed_idxs[d])

    return selectdim(frame, d, complement_idxs)
end

trimview(frame::AbstractArray, idxs::String) = trimview(frame, fits_indices(idxs))


"""
    crop(frame, shape; force_equal = true)

Crop `frame` to the size specified by `shape`, anchored by the frame center.

This will remove rows/cols of the `frame` equally on each side. When there is
an uneven difference in sizes (e.g. size 9 → 6 can't be removed equally) the
default is to increase the output size (e.g. 6 → 7) so there is equal removal
on each side. To disable this, set `force_equal=false`, which will remove the
extra slice from the end of the axis.

If `frame` is a string, it will be loaded into [`CCDData`](@ref) first.
The HDU loaded can be specified by `hdu` which by default is 1.

# Examples
```jldoctest
julia> frame = reshape(1:25, (5, 5));

julia> crop(frame, (3, 3))
3×3 Matrix{Int64}:
 7  12  17
 8  13  18
 9  14  19

julia> crop(frame, (4, 3), force_equal = false)
4×3 Matrix{Int64}:
 6  11  16
 7  12  17
 8  13  18
 9  14  19
```

# See Also
[`cropview`](@ref)
"""
crop(frame, shape; kwargs...) = copy(cropview(frame, shape; kwargs...))


"""
    cropview(frame, shape; force_equal = true)

Return a view of `frame` cropped to `shape`, anchored by the frame center.

This function is the same as [`crop`](@ref) but returns a view of the frame
rather than a copy.

!!! note
    Because this returns a view, any modification to the output will also
    modify `frame`.

# See Also
[`crop`](@ref)
"""
function cropview(frame::AbstractArray, shape; force_equal = true)
    # testing error
    if ndims(frame) != length(shape)
        throw(DimensionMismatch("Dimension mismatch between frame and shape"))
    end
    any(s -> !isa(s, Colon) && s < 1, shape) && error("crop size $shape cant't be less than 1")

    # generating idxs for cropped frame
    idxs = map(enumerate(size(frame)), shape) do (d, s1), s2
        diff = s2 isa Colon ? 0 : s1 - s2
        lower = iseven(diff) ? diff ÷ 2 : (diff - 1) ÷ 2
        upper = if isodd(diff) && force_equal
                    @warn "dimension $d changed from $s2 to $(s2 + 1)"
                    (diff - 1) ÷ 2
                elseif isodd(diff)
                    (diff + 1) ÷ 2
                else
                    diff ÷ 2
                end
        1 + lower:s1 - upper
    end

    # returning the view
    return @view frame[idxs...]
end


"""
    combine(frames...; method = median)
    combine(frames; method = median)

Combine multiple frames using `method`.

Multiple frames can also be passed in a vector or as a generator for combining.
To pass a custom method, it must have a signature like
`method(::AbstractArray; dims)`.

If `frames` are strings, they will be loaded into [`CCDData`](@ref)s first. The
HDU indices can be specified with `hdu` as either an integer or a tuple
corresponding to each file.

Header of output file (if applicable) is specified by `header_hdu` which by
default is 1.

# Examples
```jldoctest
julia> frame = [reshape(1.0:4.0, (2, 2)) for i = 1:4];

julia> combine(frame)
2×2 Matrix{Float64}:
 1.0  3.0
 2.0  4.0

julia> combine(frame, method = sum)
2×2 Matrix{Float64}:
 4.0  12.0
 8.0  16.0

```
"""
function combine(frames::Vararg{AbstractArray{<:Number}}; method = median)
    firstframe = first(frames)
    dim = ndims(firstframe) + 1
    shape = size(firstframe)
    return reshape(method(LazyStack.lazystack(frames...), dims = dim), shape)
end

combine(frames; kwargs...) = combine(frames...; kwargs...)


"""
    subtract_dark!(frame, dark_frame; data_exposure = 1, dark_exposure = 1)

In-place version of [`subtract_dark`](@ref).

# See Also
[`subtract_dark`](@ref)
"""
function subtract_dark!(frame::AbstractArray, dark_frame::AbstractArray; data_exposure = 1, dark_exposure = 1)
    factor = data_exposure / dark_exposure
    @. frame -= (dark_frame * factor)
    return frame
end


"""
    subtract_dark(frame, dark_frame; data_exposure = 1, dark_exposure = 1)

Subtract the scaled `dark_frame` from `frame`.

The dark frame is scaled by `data_exposure / dark_exposure` before subtraction
so that images taken with different exposure times can be corrected.

If either are strings, they will be loaded into [`CCDData`](@ref) first. The HDU
loaded can be specified by `hdu` as either an integer or a tuple corresponding
to each file.

# Examples
```jldoctest
julia> frame = ones(3, 3);

julia> dark_frame = ones(3, 3);

julia> subtract_dark(frame, dark_frame)
3×3 Matrix{Float64}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0

julia> subtract_dark(frame, dark_frame, data_exposure = 1, dark_exposure = 4)
3×3 Matrix{Float64}:
 0.75  0.75  0.75
 0.75  0.75  0.75
 0.75  0.75  0.75

```

# See Also
[`subtract_dark!`](@ref)
"""
function subtract_dark(frame::AbstractArray{T}, dark_frame::AbstractArray{S}; kwargs...) where {T, S}
    V = float(promote_type(T, S))
    return subtract_dark!(V.(frame), dark_frame; kwargs...)
end


"""
    gain_correct!(frame, gain)

In-place version of [`gain_correct`](@ref).

# See Also
[`gain_correct`](@ref)
"""
function gain_correct!(frame::AbstractArray, gain::Real)
    frame .*= gain
    return frame
end


"""
    gain_correct(frame, gain)

Multiply `frame` by `gain` to convert detector units (ADU) to electrons.

Correcting for gain converts the raw analogue-to-digital units (ADU) output by
the detector electronics into physical electron counts. `gain` should be
provided in electrons per ADU.

If `frame` is a string, it will be loaded into [`CCDData`](@ref) first.
The HDU loaded can be specified by `hdu` which by default is 1.

# Examples
```jldoctest
julia> frame = fill(100.0, 3, 3);

julia> gain_correct(frame, 2.5)
3×3 Matrix{Float64}:
 250.0  250.0  250.0
 250.0  250.0  250.0
 250.0  250.0  250.0

```

# See Also
[`gain_correct!`](@ref)
"""
function gain_correct(frame::AbstractArray{T}, gain::Real) where T
    return gain_correct!(float(T).(frame), gain)
end


"""
    noise_model(frame; read_noise = 0.0, gain = 1.0)

Compute a per-pixel noise estimate (standard deviation) for `frame`.

The noise model combines shot noise (Poisson noise from the photon counts) and
read noise (a fixed noise floor from the detector electronics):

```
σ(i,j) = √(max(frame[i,j], 0) / gain + (read_noise / gain)²)
```

This is the expected 1-σ uncertainty in each pixel value when the frame is
expressed in ADU. If `frame` is already in electrons (i.e., after
[`gain_correct`](@ref)), use `gain = 1`.

# Arguments
- `frame`: the image array in ADU
- `read_noise`: detector read noise in electrons (default `0.0`)
- `gain`: detector gain in electrons per ADU (default `1.0`)

# Examples
```jldoctest
julia> frame = fill(100.0, 2, 2);

julia> noise_model(frame; read_noise = 10.0, gain = 2.0)
2×2 Matrix{Float64}:
 8.66025  8.66025
 8.66025  8.66025

```
"""
function noise_model(frame::AbstractArray{T}; read_noise = 0.0, gain = 1.0) where T
    F = float(T)
    rn = F(read_noise)
    g  = F(gain)
    return @. sqrt(max(F(frame), zero(F)) / g + (rn / g)^2)
end


"""
    cosmicray_lacosmic(frame; sigma_clip = 4.5, niter = 4, read_noise = 6.5,
                              gain = 1.0, sigfrac = 0.3, objlim = 5.0,
                              background = nothing)

Detect and replace cosmic ray hits using a Laplacian edge-detection algorithm
inspired by van Dokkum (2001) (LA Cosmic).

The algorithm iterates `niter` times. In each iteration it:
1. Computes a fine-structure image using a Laplacian edge detector.
2. Estimates the noise using [`noise_model`](@ref).
3. Identifies pixels where the Laplacian-to-noise ratio exceeds `sigma_clip`
   and the contrast with the local median exceeds `objlim`.
4. Replaces flagged pixels with the local 3×3 median.

Pixels identified as cosmic ray hits are replaced by the median of their 3×3
neighbourhood. A `Bool` mask of identified cosmic rays is stored in the second
return value.

# Arguments
- `frame`: 2-D image array in ADU
- `sigma_clip`: Laplacian signal-to-noise threshold for flagging (default `4.5`)
- `niter`: number of detection–replacement iterations (default `4`)
- `read_noise`: detector read noise in electrons (default `6.5`)
- `gain`: detector gain in electrons per ADU (default `1.0`)
- `sigfrac`: fraction of `sigma_clip` used as a secondary threshold for
  neighbouring pixels (default `0.3`)
- `objlim`: minimum contrast with the fine-structure image required to flag a
  candidate pixel (default `5.0`), used to avoid flagging bright compact
  sources
- `background`: background level to subtract before noise estimation; if
  `nothing` (default) the frame median is used

# Returns
A 2-element tuple `(cleaned_frame, cosmic_ray_mask)` where `cleaned_frame`
has the same element type as `frame` and `cosmic_ray_mask` is a `BitMatrix`
that is `true` at each detected cosmic ray pixel.

# Examples
```julia
julia> frame = fill(100.0, 50, 50);

julia> frame[25, 25] = 10_000.0;  # single cosmic ray hit

julia> cleaned, mask = cosmicray_lacosmic(frame; gain = 1.5, read_noise = 5.0);

julia> mask[25, 25]
true
```
"""
function cosmicray_lacosmic(frame::AbstractArray{T};
                            sigma_clip  = 4.5,
                            niter       = 4,
                            read_noise  = 6.5,
                            gain        = 1.0,
                            sigfrac     = 0.3,
                            objlim      = 5.0,
                            background  = nothing) where T
    F = float(T)
    img = F.(frame)
    nrows, ncols = size(img)
    mask = falses(nrows, ncols)

    bg = background === nothing ? median(img) : F(background)

    for _ in 1:niter
        # 3×3 median image (used as a cleaner estimate of the true signal)
        med3 = _median3x3(img)

        # Laplacian kernel: highlights edges and point sources
        lap = _laplacian(img)

        # Noise model on the median-smoothed image
        σ = noise_model(max.(med3 .- bg, zero(F)); read_noise = read_noise, gain = gain)
        σ .= max.(σ, eps(F))

        # Fine-structure image: high Laplacian / σ ratio indicates CR
        snr = lap ./ (2 .* σ)
        snr .= max.(snr, zero(F))

        # Contrast with local median
        contrast = med3 .> zero(F)
        med5 = _median5x5(img)
        fstruct = ifelse.(contrast, img ./ max.(med3, eps(F)), zero(F))

        # Flag pixels
        new_cr = (snr .> sigma_clip) .& (fstruct .> objlim)

        # Also flag neighbours above the secondary threshold
        new_cr .|= _dilate(new_cr) .& (snr .> sigma_clip * sigfrac)

        mask .|= new_cr

        # Replace flagged pixels with 3×3 median
        for idx in findall(new_cr)
            img[idx] = med3[idx]
        end
    end

    return (convert(Array{F}, img), mask)
end

# ---- helpers for cosmicray_lacosmic ----------------------------------------

"""
    CCDReduction._median3x3(img)

Return a new array where each element is the median of the 3×3 neighbourhood
of the corresponding element in `img`. Border pixels use clamp-to-edge padding.
"""
function _median3x3(img::AbstractArray{T}) where T
    nrows, ncols = size(img)
    out = similar(img)
    buf = Vector{T}(undef, 9)
    @inbounds for j in 1:ncols, i in 1:nrows
        k = 0
        for dj in -1:1, di in -1:1
            ii = clamp(i + di, 1, nrows)
            jj = clamp(j + dj, 1, ncols)
            buf[k += 1] = img[ii, jj]
        end
        out[i, j] = median!(buf)
    end
    return out
end

"""
    CCDReduction._median5x5(img)

Return a new array where each element is the median of the 5×5 neighbourhood
of the corresponding element in `img`. Border pixels use clamp-to-edge padding.
"""
function _median5x5(img::AbstractArray{T}) where T
    nrows, ncols = size(img)
    out = similar(img)
    buf = Vector{T}(undef, 25)
    @inbounds for j in 1:ncols, i in 1:nrows
        k = 0
        for dj in -2:2, di in -2:2
            ii = clamp(i + di, 1, nrows)
            jj = clamp(j + dj, 1, ncols)
            buf[k += 1] = img[ii, jj]
        end
        out[i, j] = median!(buf)
    end
    return out
end

"""
    CCDReduction._laplacian(img)

Apply a discrete Laplacian filter to `img` and return the result, clamped to
non-negative values. Border pixels use clamp-to-edge padding.

The Laplacian kernel is:
```
 0 -1  0
-1  4 -1
 0 -1  0
```
"""
function _laplacian(img::AbstractArray{T}) where T
    nrows, ncols = size(img)
    out = similar(img)
    @inbounds for j in 1:ncols, i in 1:nrows
        up    = img[max(i-1,1),       j]
        down  = img[min(i+1,nrows),   j]
        left  = img[i,       max(j-1,1)]
        right = img[i, min(j+1,ncols)]
        center = img[i, j]
        out[i, j] = max(4*center - up - down - left - right, zero(T))
    end
    return out
end

"""
    CCDReduction._dilate(mask)

Morphological dilation of a boolean `mask` with a 3×3 structuring element.
Returns a new `BitMatrix`.
"""
function _dilate(mask::AbstractArray{Bool})
    nrows, ncols = size(mask)
    out = falses(nrows, ncols)
    @inbounds for j in 1:ncols, i in 1:nrows
        for dj in -1:1, di in -1:1
            if mask[clamp(i+di,1,nrows), clamp(j+dj,1,ncols)]
                out[i, j] = true
                break
            end
        end
    end
    return out
end
