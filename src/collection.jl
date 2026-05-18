#---------------------------------------------------------------------------------------
# ImageCollection
#---------------------------------------------------------------------------------------

"""
    CollectionRow

A single row from an [`ImageCollection`](@ref), providing named access to all
columns via `getproperty`.

The fixed columns `path`, `name`, and `hdu` are always present. Any FITS header
keywords collected by [`fitscollection`](@ref) are accessible as additional
properties.
"""
struct CollectionRow
    path :: String
    name :: String
    hdu  :: Int
    _keys :: Vector{Symbol}
    _vals :: Vector{Any}
end

function Base.getproperty(r::CollectionRow, s::Symbol)
    s === :path && return getfield(r, :path)
    s === :name && return getfield(r, :name)
    s === :hdu  && return getfield(r, :hdu)
    ks = getfield(r, :_keys)
    idx = findfirst(==(s), ks)
    idx === nothing && throw(ArgumentError("CollectionRow has no field :$s"))
    return getfield(r, :_vals)[idx]
end

function Base.show(io::IO, r::CollectionRow)
    print(io, "CollectionRow(path=$(repr(r.path)), name=$(repr(r.name)), hdu=$(r.hdu)")
    for (k, v) in zip(getfield(r, :_keys), getfield(r, :_vals))
        print(io, ", $(k)=$(repr(v))")
    end
    print(io, ")")
end

function Base.:(==)(a::CollectionRow, b::CollectionRow)
    a.path == b.path && a.name == b.name && a.hdu == b.hdu &&
    getfield(a, :_keys) == getfield(b, :_keys) &&
    getfield(a, :_vals) == getfield(b, :_vals)
end

"""
    ImageCollection

A lightweight collection of image file metadata produced by
[`fitscollection`](@ref).

An `ImageCollection` contains the file path, filename, HDU index, and any
FITS header keys found in the scanned image files. It can be iterated (each
element is a [`CollectionRow`](@ref)), indexed by integer position, filtered
with a boolean mask, and column-accessed via `getproperty`.

# Properties (column access)
- `col.paths`  — `Vector{String}` of file paths
- `col.names`  — `Vector{String}` of file names
- `col.hdus`   — `Vector{Int}` of HDU indices
- `col.KEYWORD`— `Vector{Any}` of values for FITS header keyword `KEYWORD`

# Examples
```julia
col = fitscollection("data/")

# iterate
for row in col
    println(row.path, " exptime=", row.EXPTIME)
end

# boolean filtering
science = col[col.IMAGETYP .== "LIGHT"]
```
"""
struct ImageCollection
    paths  :: Vector{String}
    names  :: Vector{String}
    hdus   :: Vector{Int}
    # header columns stored in insertion order
    _colnames :: Vector{Symbol}
    _coldata  :: Vector{Vector{Any}}
end

function ImageCollection()
    ImageCollection(String[], String[], Int[], Symbol[], Vector{Any}[])
end

Base.length(c::ImageCollection) = length(c.paths)
Base.isempty(c::ImageCollection) = isempty(c.paths)

function Base.getproperty(c::ImageCollection, s::Symbol)
    s === :paths     && return getfield(c, :paths)
    s === :names     && return getfield(c, :names)
    s === :hdus      && return getfield(c, :hdus)
    s === :_colnames && return getfield(c, :_colnames)
    s === :_coldata  && return getfield(c, :_coldata)
    cnames = getfield(c, :_colnames)
    idx = findfirst(==(s), cnames)
    idx === nothing && throw(ArgumentError("ImageCollection has no column :$s"))
    return getfield(c, :_coldata)[idx]
end

function Base.propertynames(c::ImageCollection, private::Bool = false)
    fixed = (:paths, :names, :hdus)
    extra = tuple(getfield(c, :_colnames)...)
    private ? (fixed..., extra..., :_colnames, :_coldata) : (fixed..., extra...)
end

"""
    size(collection) -> (nrows, ncols)

Return the dimensions of `collection` as `(nrows, ncols)`, where `ncols`
includes the three fixed columns (`path`, `name`, `hdu`) plus any header
keyword columns.
"""
Base.size(c::ImageCollection) = (length(c), 3 + length(getfield(c, :_colnames)))

function Base.show(io::IO, ::MIME"text/plain", c::ImageCollection)
    nr, nc = size(c)
    println(io, "ImageCollection with $nr rows and $nc columns")
    print(io, "  Fixed columns : path, name, hdu")
    cnames = getfield(c, :_colnames)
    if !isempty(cnames)
        print(io, "\n  Header columns: ", join(cnames, ", "))
    end
end

Base.show(io::IO, c::ImageCollection) = show(io, MIME"text/plain"(), c)

"""
    getindex(collection, i)

Return the `i`-th row of `collection` as a [`CollectionRow`](@ref).
"""
function Base.getindex(c::ImageCollection, i::Integer)
    cnames = getfield(c, :_colnames)
    cdata  = getfield(c, :_coldata)
    vals = [cdata[j][i] for j in eachindex(cnames)]
    return CollectionRow(c.paths[i], c.names[i], c.hdus[i], copy(cnames), vals)
end

"""
    getindex(collection, mask::AbstractVector{Bool})

Return a new [`ImageCollection`](@ref) with only the rows where `mask` is `true`.
"""
function Base.getindex(c::ImageCollection, mask::AbstractVector{Bool})
    length(mask) == length(c) || throw(DimensionMismatch(
        "mask length ($(length(mask))) does not match collection length ($(length(c)))"))
    cnames = getfield(c, :_colnames)
    cdata  = getfield(c, :_coldata)
    ImageCollection(
        c.paths[mask],
        c.names[mask],
        c.hdus[mask],
        copy(cnames),
        [col[mask] for col in cdata],
    )
end

"""
    iterate(collection)

Iterate over the rows of `collection`, yielding [`CollectionRow`](@ref)
objects.
"""
function Base.iterate(c::ImageCollection, state = 1)
    state > length(c) && return nothing
    return (c[state], state + 1)
end

Base.eltype(::Type{ImageCollection}) = CollectionRow

#---------------------------------------------------------------------------------------
# Utility helpers
#---------------------------------------------------------------------------------------

# parses the name and returns it with or without extension
parse_name(filename, ext::AbstractString, ::Val{false}) = first(rsplit(filename, ext, limit=2))

function parse_name(filename, ext::Regex, ::Val{false})
    idxs = findall(ext, filename)
    return filename[1:first(last(idxs)) - 1]
end

parse_name(filename, ext, ::Val{true}) = filename

# utility function for generating filename
function generate_filename(path, save_location, save_prefix, save_suffix, save_delim, ext)
    filename = basename(path)
    modified_name, extension = parse_name_ext(filename, "." * ext)
    if !isnothing(save_prefix)
        modified_name = string(save_prefix, save_delim, modified_name)
    end
    if !isnothing(save_suffix)
        modified_name = string(modified_name, save_delim, save_suffix)
    end
    file_path = joinpath(save_location, modified_name * extension)
    return file_path
end

# utility function to return filename and extension separately
# returns extension including "." at the beginning
function parse_name_ext(filename, ext)
    idxs = findall(ext, filename)
    length(idxs) == 0 && return (filename, "")
    breaking_index = first(last(idxs))
    return filename[1:breaking_index - 1], filename[breaking_index:end]
end

#---------------------------------------------------------------------------------------
# fitscollection
#---------------------------------------------------------------------------------------

@doc raw"""
    fitscollection(dir;
                   recursive    = true,
                   abspath      = true,
                   keepext      = true,
                   ext          = r"fits(\.tar\.gz)?",
                   exclude      = nothing,
                   exclude_dir  = nothing,
                   exclude_key  = ("", "HISTORY"))

Walk through `dir` collecting FITS files, scanning their headers, and
returning an [`ImageCollection`](@ref) that can be used with the iterators
[`arrays`](@ref), [`filenames`](@ref), and [`ccds`](@ref) for batch processing.

If `recursive` is `false`, no subdirectories will be walked. The collection
contains the fixed columns `path`, `name`, and `hdu`, plus one column per
distinct FITS header keyword found across all scanned files.  Missing values
(a keyword present in some but not all files) are stored as `nothing`.

!!! note "Duplicate Keys"
    In certain cases there are multiple FITS headers with the same key, e.g.,
    `COMMENT`. In these cases only the first instance of the key-value pair
    will be stored.

If `abspath` is `true`, the path in the collection will be absolute. If
`keepext` is `true`, the name will include the file extension, given by `ext`.
`ext` is used with `endswith` to filter for FITS files compatible with
`FITSIO.FITS`.

`exclude` is a pattern used with `occursin` to exclude certain filenames.
Similarly, `exclude_dir` allows excluding entire folders.
`exclude_key` lists header keywords to skip when building the collection
(e.g. `("", "HISTORY")` skips empty keys and `HISTORY` cards).

# Examples
```julia
col = fitscollection("data/")

# access columns
col.paths
col.EXPTIME  # all EXPTIME values as a Vector

# filter to science frames only
sci = col[col.IMAGETYP .== "LIGHT"]

# iterate
for row in col
    img = CCDData(row.path; hdu = row.hdu)
end
```
"""
function fitscollection(basedir::String;
                        recursive   = true,
                        abspath     = true,
                        keepext     = true,
                        ext         = r"fits(\.tar\.gz)?"i,
                        exclude     = nothing,
                        exclude_dir = nothing,
                        exclude_key = ("", "HISTORY"))

    # Accumulate rows as a vector of NamedTuples temporarily
    # then build ImageCollection columns
    all_paths  = String[]
    all_names  = String[]
    all_hdus   = Int[]

    # Each row's header data stored as Dict{Symbol,Any}
    row_headers = Dict{Symbol,Any}[]

    # Track all header keys seen (in order of first appearance)
    header_key_order = Symbol[]
    header_key_set   = Set{Symbol}()

    for (root, dirs, files) in walkdir(basedir)
        recursive || root == basedir || continue
        if exclude_dir !== nothing
            occursin(exclude_dir, root) && continue
        end
        for filename in files
            endswith(filename, ext) || continue
            if exclude !== nothing
                occursin(exclude, filename) && continue
            end
            location = joinpath(root, filename)
            fits_data = FITS(location)

            for (index, hdu) in enumerate(fits_data)
                hdu isa ImageHDU || continue
                header_data = read_header(hdu)
                path = abspath ? Base.abspath(location) : location
                name = parse_name(filename, "." * ext, Val(keepext))

                _keys  = filter(k -> k ∉ exclude_key, keys(header_data))
                unique_inds = unique(idx -> _keys[idx], eachindex(_keys))
                unique_keys = _keys[unique_inds]

                push!(all_paths, path)
                push!(all_names, name)
                push!(all_hdus,  index)

                row_hdr = Dict{Symbol,Any}()
                for k in unique_keys
                    sk = Symbol(k)
                    row_hdr[sk] = header_data[k]
                    if sk ∉ header_key_set
                        push!(header_key_set, sk)
                        push!(header_key_order, sk)
                    end
                end
                push!(row_headers, row_hdr)
            end
            close(fits_data)
        end
    end

    # Build column vectors (fill missing with nothing)
    n = length(all_paths)
    coldata = [Vector{Any}(nothing, n) for _ in header_key_order]

    for (i, row_hdr) in enumerate(row_headers)
        for (ci, k) in enumerate(header_key_order)
            if haskey(row_hdr, k)
                coldata[ci][i] = row_hdr[k]
            end
        end
    end

    return ImageCollection(all_paths, all_names, all_hdus,
                           header_key_order, coldata)
end

#---------------------------------------------------------------------------------------
# Iterators over collections
#---------------------------------------------------------------------------------------

"""
    arrays(collection)

Return a lazy iterator that loads each image in `collection` as an `Array`.

Each call to `next` opens the FITS file at `row.path`, reads the HDU at
`row.hdu`, and closes the file handle.

# Examples
```julia
col = fitscollection("data/")
for arr in arrays(col)
    println(size(arr))
end
```
"""
function arrays(collection::ImageCollection)
    return (begin
        fh  = FITS(row.path)
        arr = getdata(fh[row.hdu])
        close(fh)
        arr
    end for row in collection)
end


"""
    filenames(collection)

Return a lazy iterator over the file paths in `collection`.

# Examples
```julia
col = fitscollection("data/")
for path in filenames(col)
    println(path)
end
```
"""
filenames(collection::ImageCollection) = (row.path for row in collection)


"""
    ccds(collection)

Return a lazy iterator that loads each image in `collection` as a
[`CCDData`](@ref).

# Examples
```julia
col = fitscollection("data/")
for ccd in ccds(col)
    println(size(ccd))
end
```
"""
ccds(collection::ImageCollection) = (CCDData(row.path; hdu = row.hdu) for row in collection)


#---------------------------------------------------------------------------------------
# f-form iterators with optional saving
#---------------------------------------------------------------------------------------

"""
    ccds(f, collection;
         path         = nothing,
         save_prefix  = nothing,
         save_suffix  = nothing,
         save         = any(!isnothing, (save_prefix, path, save_suffix)),
         save_delim   = "_",
         ext          = r"fits(\\.tar\\.gz)?"i)

Iterate over the [`CCDData`](@ref) frames in `collection`, apply `f` to each,
and return a `Vector` of results.

The outputs from `f` can optionally be saved as FITS files.  `save_prefix`
adds a prefix to each filename separated by `save_delim`; `save_suffix` adds a
suffix before the extension.  Files are saved in their original directory
unless `path` is given.  Saving is activated automatically if any of
`save_prefix`, `path`, or `save_suffix` is set, but can be overridden with
`save`.

# Example
```julia
col = fitscollection("data/")

# process and save
results = ccds(col; path = "output/", save_prefix = "reduced") do img
    trim(img, (:, 513:524))
end
```
"""
function ccds(f,
              collection::ImageCollection;
              path         = nothing,
              save_prefix  = nothing,
              save_suffix  = nothing,
              save         = any(!isnothing, (save_prefix, path, save_suffix)),
              save_delim   = "_",
              ext          = r"fits(\.tar\.gz)?"i)

    map(collection) do row
        img = CCDData(row.path; hdu = row.hdu)
        out = f(img)
        if save
            save_path = generate_filename(row.path,
                                          something(path, dirname(row.path)),
                                          save_prefix, save_suffix, save_delim, ext)
            writefits(save_path, out)
        end
        out
    end
end


"""
    filenames(f, collection;
              path         = nothing,
              save_prefix  = nothing,
              save_suffix  = nothing,
              save         = any(!isnothing, (save_prefix, path, save_suffix)),
              save_delim   = "_",
              ext          = r"fits(\\.tar\\.gz)?"i)

Iterate over the file paths in `collection`, apply `f` to each path, and
return a `Vector` of results.  Saving behaviour is identical to the `f`-form
of [`ccds`](@ref).

# Example
```julia
col = fitscollection("data/")
sizes = filenames(col) do path
    FITS(f -> size(read(f[1])), path)
end
```
"""
function filenames(f,
                   collection::ImageCollection;
                   path         = nothing,
                   save_prefix  = nothing,
                   save_suffix  = nothing,
                   save         = any(!isnothing, (save_prefix, path, save_suffix)),
                   save_delim   = "_",
                   ext          = r"fits(\.tar\.gz)?"i)

    map(collection) do row
        out = f(row.path)
        if save
            save_path = generate_filename(row.path,
                                          something(path, dirname(row.path)),
                                          save_prefix, save_suffix, save_delim, ext)
            writefits(save_path, out)
        end
        out
    end
end


"""
    arrays(f, collection;
           path         = nothing,
           save_prefix  = nothing,
           save_suffix  = nothing,
           save         = any(!isnothing, (save_prefix, path, save_suffix)),
           save_delim   = "_",
           ext          = r"fits(\\.tar\\.gz)?"i)

Iterate over the image arrays in `collection`, apply `f` to each, and return
a `Vector` of results.  Saving behaviour is identical to the `f`-form of
[`ccds`](@ref).

# Example
```julia
col = fitscollection("data/")
trimmed = arrays(col; path = "output/", save_prefix = "trim") do arr
    trim(arr, (:, 513:524))
end
```
"""
function arrays(f,
                collection::ImageCollection;
                path         = nothing,
                save_prefix  = nothing,
                save_suffix  = nothing,
                save         = any(!isnothing, (save_prefix, path, save_suffix)),
                save_delim   = "_",
                ext          = r"fits(\.tar\.gz)?"i)

    map(collection) do row
        fh  = FITS(row.path)
        arr = getdata(fh[row.hdu])
        close(fh)
        out = f(arr)
        if save
            save_path = generate_filename(row.path,
                                          something(path, dirname(row.path)),
                                          save_prefix, save_suffix, save_delim, ext)
            writefits(save_path, out)
        end
        out
    end
end
