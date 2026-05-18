using CCDReduction: parse_name, generate_filename, parse_name_ext, writefits

# ---------------------------------------------------------------------------
# ImageCollection / fitscollection
# ---------------------------------------------------------------------------
@testset "fitscollection – basic" begin
    col = fitscollection(TEST_DATA_DIR)

    # Each of the 5 mock FITS files should produce one row
    @test length(col) == 5
    nr, nc = size(col)
    @test nr == 5
    @test nc >= 3  # at least path, name, hdu + header columns

    # Fixed columns accessible as properties
    @test col.paths isa Vector{String}
    @test col.names isa Vector{String}
    @test col.hdus  isa Vector{Int}
    @test all(col.hdus .== 1)

    # Header keyword columns (EXPTIME, IMAGETYP, FILTER were written for some)
    @test :EXPTIME   ∈ propertynames(col)
    @test :IMAGETYP  ∈ propertynames(col)

    # names include extension by default
    @test all(endswith.(col.names, ".fits"))

    # keepext = false
    col2 = fitscollection(TEST_DATA_DIR; keepext = false)
    @test all(!endswith(n, ".fits") for n in col2.names)
end

@testset "fitscollection – exclude options" begin
    # Use a fresh temp dir to avoid polluting TEST_DATA_DIR
    mktempdir() do tmpdir
        cp(FRAME_A_PATH, joinpath(tmpdir, "frameA.fits"))
        cp(FRAME_B_PATH, joinpath(tmpdir, "frameB.fits"))
        cp(BIAS_PATH,    joinpath(tmpdir, "bias.fits"))

        # exclude by filename pattern
        col_ex = fitscollection(tmpdir; exclude = "bias")
        @test all(!occursin("bias", n) for n in col_ex.names)
        @test length(col_ex) == 2

        # exclude_dir
        col_no_dir = fitscollection(tmpdir; exclude_dir = tmpdir)
        @test length(col_no_dir) == 0

        # non-recursive: subdir file should not appear without recursive=true
        subdir = mktempdir(tmpdir)
        make_fits(joinpath(subdir, "sub.fits"), fill(1.0, 4, 4))
        col_rec   = fitscollection(tmpdir; recursive = true)
        col_norec = fitscollection(tmpdir; recursive = false)
        @test length(col_rec) == length(col_norec) + 1
    end
end

@testset "ImageCollection – indexing & iteration" begin
    col = fitscollection(TEST_DATA_DIR)

    # Integer indexing returns CollectionRow
    row = col[1]
    @test row isa CCDReduction.CollectionRow
    @test row.path isa String
    @test row.hdu == 1
    @test row.name isa String

    # Property access on row
    @test row.EXPTIME isa Number
    @test row.IMAGETYP isa String

    # Boolean mask indexing
    mask = col.IMAGETYP .== "LIGHT"
    col_light = col[mask]
    @test col_light isa ImageCollection
    @test length(col_light) == 2   # frameA and frameB are LIGHT frames

    # Iteration
    paths_iter = [row.path for row in col]
    @test sort(paths_iter) == sort(col.paths)

    # size
    @test size(col)[1] == length(col)
end

# ---------------------------------------------------------------------------
# arrays / filenames / ccds iterators
# ---------------------------------------------------------------------------
@testset "arrays iterator" begin
    col = fitscollection(TEST_DATA_DIR)
    arrs = collect(arrays(col))
    @test length(arrs) == length(col)
    @test all(a -> a isa Array, arrs)

    # Verify values round-trip for frameA specifically
    mktempdir() do tmpdir
        cp(FRAME_A_PATH, joinpath(tmpdir, "frameA.fits"))
        col_a = fitscollection(tmpdir)
        for arr in arrays(col_a)
            @test size(arr) == size(FRAME_A_DATA)
            @test arr ≈ FRAME_A_DATA
        end
    end
end

@testset "filenames iterator" begin
    col  = fitscollection(TEST_DATA_DIR)
    fns  = collect(filenames(col))
    @test fns == col.paths
end

@testset "ccds iterator" begin
    col = fitscollection(TEST_DATA_DIR)
    for ccd in ccds(col)
        @test ccd isa CCDData
    end
end

# ---------------------------------------------------------------------------
# f-form iterators with saving
# ---------------------------------------------------------------------------
@testset "ccds(f, col) – saving" begin
    # Use an isolated dir with only frameA
    mktempdir() do srcdir
        cp(FRAME_A_PATH, joinpath(srcdir, "frameA.fits"))
        col = fitscollection(srcdir)

        mktempdir() do savedir
            results = ccds(col; path = savedir, save_prefix = "reduced", save_suffix = "v1") do img
                trim(img, (:, 60:64))
            end

            @test length(results) == 1
            @test results[1] isa CCDData
            @test size(results[1]) == (64, 59)

            col2 = fitscollection(savedir; recursive = false)
            @test length(col2) == 1
            @test occursin("reduced", col2[1].name)
            @test occursin("v1",      col2[1].name)

            # Data round-trip
            arr = first(arrays(col2))
            @test arr ≈ results[1].data
        end
    end
end

@testset "arrays(f, col) – saving" begin
    mktempdir() do srcdir
        cp(FRAME_A_PATH, joinpath(srcdir, "frameA.fits"))
        col = fitscollection(srcdir)

        mktempdir() do savedir
            results = arrays(col; path = savedir, save_prefix = "arr_out") do arr
                trim(arr, (:, 60:64))
            end

            @test results[1] ≈ trim(FRAME_A_DATA, (:, 60:64))

            col2 = fitscollection(savedir; recursive = false)
            arr2 = first(arrays(col2))
            @test arr2 ≈ results[1]
        end
    end
end

@testset "filenames(f, col) – saving" begin
    mktempdir() do srcdir
        cp(FRAME_A_PATH, joinpath(srcdir, "frameA.fits"))
        col = fitscollection(srcdir)

        mktempdir() do savedir
            results = filenames(col; path = savedir, save_prefix = "fn_out") do path
                arr = getdata(FITS(path)[1])
                arr .+ 1.0f0
            end

            @test results[1] ≈ FRAME_A_DATA .+ 1.0f0

            col2 = fitscollection(savedir; recursive = false)
            arr2 = first(arrays(col2))
            @test arr2 ≈ results[1]
        end
    end
end

@testset "save without explicit path (same dir as input)" begin
    # copy one file to a fresh temp dir so we don't pollute TEST_DATA_DIR
    tmpdir = mktempdir()
    cp(FRAME_A_PATH, joinpath(tmpdir, "frameA.fits"))
    col = fitscollection(tmpdir)
    @test length(col) == 1

    _ = ccds(col; save_prefix = "inplace") do img; img; end
    col2 = fitscollection(tmpdir)
    @test length(col2) == 2
end

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------
@testset "parse_name" begin
    @test parse_name("abc.fits", "."*"fits",            Val(true))  == "abc.fits"
    @test parse_name("abc.fits.tar.gz", "."*"fits.tar.gz", Val(false)) == "abc"
    @test parse_name("foo.fits.fits",   "."*"fits",         Val(false)) == "foo.fits"
    @test parse_name("foo.fits",        "."*r"fits(\.tar\.gz)?"i, Val(false)) == "foo"
end

@testset "generate_filename" begin
    @test generate_filename("home/abcd.fits", @__DIR__, "t1", "t2", "_",  r"fits(\.tar\.gz)?"i) == joinpath(@__DIR__, "t1_abcd_t2.fits")
    @test generate_filename("home/abcd.fits", @__DIR__, nothing, "t2", "_",  r"fits(\.tar\.gz)?"i) == joinpath(@__DIR__, "abcd_t2.fits")
    @test generate_filename("home/abcd.fits", @__DIR__, "t1", nothing, "_",  r"fits(\.tar\.gz)?"i) == joinpath(@__DIR__, "t1_abcd.fits")
    @test generate_filename("home/abcd.fits", @__DIR__, nothing, nothing, "_", r"fits(\.tar\.gz)?"i) == joinpath(@__DIR__, "abcd.fits")
    @test generate_filename("home/abcd.fits", @__DIR__, "t1", nothing, "__", r"fits(\.tar\.gz)?"i) == joinpath(@__DIR__, "t1__abcd.fits")
end

@testset "writefits round-trip" begin
    mktempdir() do dir
        path   = joinpath(dir, "test.fits")
        data   = rand(5, 10)
        writefits(path, data)
        fh     = FITS(path)
        result = getdata(fh[1])
        @test result ≈ data
        close(fh)
    end
end
