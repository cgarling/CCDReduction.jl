using CCDReduction: getdata, writefits

function test_header(ccd1::CCDData, ccd2::CCDData)
    h1 = ccd1.hdr
    h2 = ccd2.hdr
    @test keys(h1) == keys(h2)
    for (k1, k2) in zip(keys(h1), keys(h2))
        @test h1[k1] == h2[k2]
    end
end

# ---------------------------------------------------------------------------
# getdata / writefits round-trip
# ---------------------------------------------------------------------------
@testset "getdata / writefits" begin
    mktempdir() do dir
        path = joinpath(dir, "rw.fits")
        data = rand(12, 15)
        writefits(path, data)
        fh = FITS(path)
        @test getdata(fh[1]) ≈ data
        close(fh)

        # CCDData round-trip
        path2 = joinpath(dir, "rw2.fits")
        ccd = CCDData(fill(3.0, 8, 8))
        writefits(path2, ccd)
        fh2 = FITS(path2)
        @test getdata(fh2[1]) ≈ ccd.data
        close(fh2)
    end
end

# ---------------------------------------------------------------------------
# default_header
# ---------------------------------------------------------------------------
@testset "default_header / CCDData constructors" begin
    # CCDData from raw array should get a valid default header
    ccd = CCDData(zeros(4, 4))
    @test ccd isa CCDData
    @test ccd[:SIMPLE] == true

    # CCDData from FITS path
    ccd2 = CCDData(FRAME_A_PATH; hdu = 1)
    @test ccd2 isa CCDData
    @test size(ccd2) == (64, 64)

    # CCDData from ImageHDU
    FITS(FRAME_A_PATH) do fh
        ccd3 = CCDData(fh[1])
        @test size(ccd3) == (64, 64)
    end
end

# ---------------------------------------------------------------------------
# subtract_bias with FITS inputs
# ---------------------------------------------------------------------------
@testset "subtract_bias (FITS)" begin
    hdu_frame = CCDData(FRAME_A_PATH)
    hdu_bias  = CCDData(BIAS_PATH)
    arr_frame = getdata(FITS(FRAME_A_PATH)[1])
    arr_bias  = getdata(FITS(BIAS_PATH)[1])

    # CCDData - CCDData
    result = subtract_bias(hdu_frame, hdu_bias)
    @test result isa CCDData
    @test result.data ≈ arr_frame .- arr_bias
    test_header(result, hdu_frame)

    # Array - CCDData
    result2 = subtract_bias(arr_frame, hdu_bias)
    @test result2 isa Array
    @test result2 ≈ arr_frame .- arr_bias

    # String - CCDData
    result3 = subtract_bias(FRAME_A_PATH, hdu_bias; hdu = 1)
    @test result3 isa CCDData
    @test result3.data ≈ arr_frame .- arr_bias

    # CCDData - String
    result4 = subtract_bias(hdu_frame, BIAS_PATH; hdu = 1)
    @test result4 isa CCDData
    @test result4.data ≈ arr_frame .- arr_bias

    # Mutating: CCDData - CCDData
    hdu_frame2 = CCDData(FRAME_A_PATH)
    subtract_bias!(hdu_frame2, hdu_bias)
    @test hdu_frame2.data ≈ arr_frame .- arr_bias

    # Integer frame: InexactError when bias is float with fractional part
    hdu_int = CCDData(fill(Int32(2), 5, 5))
    hdu_frac_bias = CCDData(fill(2.5, 5, 5))
    @test_throws InexactError subtract_bias!(hdu_int, hdu_frac_bias)
end

# ---------------------------------------------------------------------------
# subtract_overscan with FITS inputs
# ---------------------------------------------------------------------------
@testset "subtract_overscan (FITS)" begin
    hdu_frame = CCDData(FRAME_A_PATH)
    arr_frame = Float64.(getdata(FITS(FRAME_A_PATH)[1]))

    # CCDData with tuple indices
    result = subtract_overscan(hdu_frame, (:, 60:64))
    @test result isa CCDData
    @test result.data ≈ subtract_overscan(arr_frame, (:, 60:64))

    # String
    result2 = subtract_overscan(FRAME_A_PATH, (:, 60:64); hdu = 1)
    @test result2 isa CCDData
    @test result2.data ≈ subtract_overscan(arr_frame, (:, 60:64))

    # Mutating
    hdu_frame2 = CCDData(FRAME_A_PATH)
    subtract_overscan!(hdu_frame2, (:, 60:64))
    @test hdu_frame2.data ≈ subtract_overscan(arr_frame, (:, 60:64))
end

# ---------------------------------------------------------------------------
# flat_correct with FITS inputs
# ---------------------------------------------------------------------------
@testset "flat_correct (FITS)" begin
    hdu_frame = CCDData(FRAME_A_PATH)
    hdu_flat  = CCDData(FLAT_PATH)
    arr_frame = Float64.(getdata(FITS(FRAME_A_PATH)[1]))
    arr_flat  = Float64.(getdata(FITS(FLAT_PATH)[1]))

    # CCDData - CCDData (flat = 1.0, so result = frame)
    result = flat_correct(hdu_frame, hdu_flat)
    @test result isa CCDData
    @test result.data ≈ arr_frame   # flat is all-ones, mean=1, so no change

    # With norm_value = 1
    result2 = flat_correct(hdu_frame, hdu_flat; norm_value = 1.0)
    @test result2.data ≈ arr_frame

    # String - CCDData
    result3 = flat_correct(FRAME_A_PATH, hdu_flat; norm_value = 1.0, hdu = 1)
    @test result3 isa CCDData
    @test result3.data ≈ arr_frame

    # Mutating
    hdu_frame2 = CCDData(copy(Float64.(hdu_frame.data)), hdu_frame.hdr)
    flat_correct!(hdu_frame2, hdu_flat; norm_value = 1.0)
    @test hdu_frame2.data ≈ arr_frame

    # Type promotion: integer frame with float flat
    hdu_int_frame = CCDData(fill(Int32(4), 5, 5))
    hdu_float_flat = CCDData(fill(2.0, 5, 5))
    result4 = flat_correct(hdu_int_frame, hdu_float_flat; norm_value = 1.0)
    @test result4 isa CCDData
    @test result4.data ≈ fill(2.0, 5, 5)
end

# ---------------------------------------------------------------------------
# trim with FITS inputs
# ---------------------------------------------------------------------------
@testset "trim (FITS)" begin
    hdu_frame = CCDData(FRAME_A_PATH)
    arr_frame = getdata(FITS(FRAME_A_PATH)[1])

    # CCDData tuple indices
    result = trim(hdu_frame, (:, 60:64))
    @test result isa CCDData
    @test result.data == trim(arr_frame, (:, 60:64))

    # String
    result2 = trim(FRAME_A_PATH, (:, 60:64); hdu = 1)
    @test result2 isa CCDData
    @test result2.data == trim(arr_frame, (:, 60:64))

    # trimview returns a view
    v = trimview(hdu_frame, (:, 60:64))
    @test v isa CCDData
    @test v.data isa SubArray
end

# ---------------------------------------------------------------------------
# crop with FITS inputs
# ---------------------------------------------------------------------------
@testset "crop (FITS)" begin
    hdu_frame = CCDData(FRAME_A_PATH)
    arr_frame = getdata(FITS(FRAME_A_PATH)[1])

    result = crop(hdu_frame, (32, 32))
    @test result isa CCDData
    @test result.data == crop(arr_frame, (32, 32))

    result2 = crop(FRAME_A_PATH, (32, 32); hdu = 1)
    @test result2 isa CCDData
    @test result2.data == crop(arr_frame, (32, 32))
end

# ---------------------------------------------------------------------------
# combine with FITS inputs
# ---------------------------------------------------------------------------
@testset "combine (FITS)" begin
    hdu = CCDData(FRAME_A_PATH)
    arr = getdata(FITS(FRAME_A_PATH)[1])

    # Vector of CCDData
    result = combine([hdu, hdu, hdu])
    @test result isa CCDData
    @test result.data ≈ combine([arr, arr, arr])

    # Varargs CCDData
    result2 = combine(hdu, hdu, hdu)
    @test result2.data ≈ result.data

    # header_hdu selection
    ccd1 = CCDData(ones(5, 5))
    ccd2 = CCDData(fill(2.0, 5, 5))
    ccd2.hdr["SIMPLE"] = false
    r = combine(ccd1, ccd2; header_hdu = 2, method = sum)
    @test r.hdr["SIMPLE"] == false

    # String inputs
    result3 = combine([FRAME_A_PATH, FRAME_A_PATH]; method = sum)
    @test result3 isa CCDData
    @test result3.data ≈ 2 .* Float64.(arr)
end

# ---------------------------------------------------------------------------
# subtract_dark with FITS inputs
# ---------------------------------------------------------------------------
@testset "subtract_dark (FITS)" begin
    hdu_frame = CCDData(FRAME_A_PATH)
    hdu_dark  = CCDData(DARK_PATH)
    arr_frame = Float64.(getdata(FITS(FRAME_A_PATH)[1]))
    arr_dark  = Float64.(getdata(FITS(DARK_PATH)[1]))

    # Equal exposures
    result = subtract_dark(hdu_frame, hdu_dark;
                           data_exposure = 1.0, dark_exposure = 1.0)
    @test result isa CCDData
    @test result.data ≈ arr_frame .- arr_dark

    # String inputs
    result2 = subtract_dark(FRAME_A_PATH, DARK_PATH;
                             data_exposure = 1.0, dark_exposure = 1.0, hdu = 1)
    @test result2 isa CCDData
    @test result2.data ≈ arr_frame .- arr_dark

    # Mutating
    hdu_frame2 = CCDData(FRAME_A_PATH)
    subtract_dark!(hdu_frame2, hdu_dark; data_exposure = 1.0, dark_exposure = 1.0)
    @test hdu_frame2.data ≈ arr_frame .- arr_dark
end

# ---------------------------------------------------------------------------
# gain_correct with FITS inputs
# ---------------------------------------------------------------------------
@testset "gain_correct (FITS)" begin
    hdu_frame = CCDData(FRAME_A_PATH)
    arr_frame = Float64.(getdata(FITS(FRAME_A_PATH)[1]))

    result = gain_correct(hdu_frame, 2.0)
    @test result isa CCDData
    @test result.data ≈ 2.0 .* arr_frame

    result2 = gain_correct(FRAME_A_PATH, 2.0; hdu = 1)
    @test result2 isa CCDData
    @test result2.data ≈ 2.0 .* arr_frame
end

# ---------------------------------------------------------------------------
# Header access via CCDData
# ---------------------------------------------------------------------------
@testset "CCDData header access" begin
    ccd = CCDData(zeros(5, 5))
    @test ccd[:SIMPLE] == true
    @test ccd[:SIMPLE] == ccd["SIMPLE"]
    ccd[:SIMPLE] = false
    @test ccd["SIMPLE"] == false
end
