using CCDReduction: axes_min_length, fits_indices, convert_value

# ---------------------------------------------------------------------------
# subtract_bias
# ---------------------------------------------------------------------------
@testset "subtract_bias – correctness" begin
    # Float: basic subtraction
    frame = [1.0 2.2 3.3 4.5]
    bias  = [0.0 0.2 0.3 0.5]
    @test subtract_bias(frame, bias) ≈ [1.0 2.0 3.0 4.0]

    # Integer frames
    frame_i = [10  20  30; 40  50  60]
    bias_i  = [ 5  10  15;  5  10  15]
    @test subtract_bias(frame_i, bias_i) == [5 10 15; 35 40 45]

    # Uniform bias (common case)
    @test subtract_bias(ones(500, 500), ones(500, 500)) == zeros(500, 500)

    # Mutating version
    f = copy(frame)
    subtract_bias!(f, bias)
    @test f ≈ [1.0 2.0 3.0 4.0]

    # Non-mutating does not modify input
    orig = fill(5.0, 4, 4)
    _ = subtract_bias(orig, fill(1.0, 4, 4))
    @test orig == fill(5.0, 4, 4)

    # Type-stability
    @test @inferred(subtract_bias(ones(4, 4), ones(4, 4))) == zeros(4, 4)
    @test @inferred(subtract_bias!(ones(4, 4), ones(4, 4))) == zeros(4, 4)

    # Error on size mismatch
    @test_throws DimensionMismatch subtract_bias(ones(5, 1), ones(5, 5))
end

# ---------------------------------------------------------------------------
# subtract_overscan
# ---------------------------------------------------------------------------
@testset "subtract_overscan – correctness" begin
    frame = Float64[4.0 2.0 3.0 1.0 1.0]

    # Tuple indices
    @test subtract_overscan(frame, (:, 4:5), dims = 2) ≈ [3.0 1.0 2.0 0.0 0.0]

    # FITS-style string
    @test subtract_overscan(frame, "[4:5, 1:1]", dims = 2) ≈ [3.0 1.0 2.0 0.0 0.0]

    # 2-D frame: overscan along rows
    f2d = ones(5, 10)
    f2d[:, 9:10] .= 0.5
    result = subtract_overscan(f2d, (:, 9:10))
    @test result[:, 9:10] ≈ zeros(5, 2)
    @test result[:, 1:8]  ≈ fill(0.5, 5, 8)

    # Integer frame (overscan should round correctly)
    f_int = Int16[10 10 10 2 2; 10 10 10 2 2]
    r_int = subtract_overscan(f_int, (:, 4:5))
    @test r_int[:, 1:3] == fill(Int16(8), 2, 3)

    # Mutating version
    fm = [1.0 1.0 1.0 0.5 0.5]
    subtract_overscan!(fm, (:, 4:5))
    @test fm[1, 1:3] ≈ [0.5, 0.5, 0.5]

    # Type stability
    @test @inferred(subtract_overscan(ones(5, 10), (:, 9:10))) isa Matrix{Float64}

    # Out-of-bounds error
    @test_throws BoundsError subtract_overscan(ones(5, 6), (5:7, :))
end

# ---------------------------------------------------------------------------
# flat_correct
# ---------------------------------------------------------------------------
@testset "flat_correct – correctness" begin
    frame = ones(3, 3)
    flat  = fill(2.0, 3, 3)

    # With norm_value = mean(flat) = 2, result should be ones
    @test flat_correct(frame, flat) ≈ ones(3, 3)

    # With custom norm_value
    @test flat_correct(frame, flat, norm_value = 1.0) ≈ fill(0.5, 3, 3)
    @test flat_correct(frame, flat, norm_value = 4.0) ≈ fill(2.0, 3, 3)

    # Integer frame promoted to float
    frame_i = fill(Int32(4), 5, 5)
    flat_i  = fill(Int32(2), 5, 5)
    result  = flat_correct(frame_i, flat_i, norm_value = 1.0)
    @test result isa Matrix{Float64}
    @test result ≈ fill(2.0, 5, 5)

    # Flat normalised by its mean: science / (flat/mean(flat)) = science when flat is uniform
    @test flat_correct(fill(3.0, 4, 4), fill(3.0, 4, 4)) ≈ fill(3.0, 4, 4)

    # Mutating version
    fm = ones(3, 3)
    flat_correct!(fm, flat; norm_value = 1.0)
    @test fm ≈ fill(0.5, 3, 3)

    # Non-mutating does not modify input
    orig = fill(10.0, 4, 4)
    _ = flat_correct(orig, fill(2.0, 4, 4))
    @test orig == fill(10.0, 4, 4)

    # Type-stability
    @test @inferred(flat_correct(ones(4, 4), ones(4, 4))) isa Matrix{Float64}

    # Error: negative norm_value
    @test_throws ErrorException flat_correct(ones(4, 4), ones(4, 4), norm_value = -1.0)
    # Error: size mismatch
    @test_throws DimensionMismatch flat_correct(ones(5, 5), ones(5, 6))
end

# ---------------------------------------------------------------------------
# trim / trimview
# ---------------------------------------------------------------------------
@testset "trim – correctness" begin
    x = reshape(1:25, 5, 5)

    # trim columns
    @test trim(x, (:, 4:5)) == [1:5 6:10 11:15]

    # trim rows
    @test trim(x, (1:2, :)) == [3:5 8:10 13:15 18:20 23:25]

    # FITS string equivalence
    @test trim(x, "[4:5, :]") == trim(x, (:, 4:5))
    @test trim(x, "[:, 4:5]") == trim(x, (4:5, :))

    # trim returns a copy; trimview returns a view
    @test trim(ones(5, 5), (:, 3:5)) isa Array
    @test trimview(ones(5, 5), (:, 3:5)) isa SubArray

    # View shares memory
    y = collect(x)
    v = trimview(y, (:, 3:5))
    v[1] = 999
    @test y[1, 1] == 999

    # Float frame
    f = fill(1.5, 6, 6)
    @test trim(f, (:, 6:6)) == fill(1.5, 6, 5)

    # Errors
    @test_throws ErrorException trim(ones(5, 5), (4:5, 1:4))  # both dims constrained
    @test_throws ErrorException trim(ones(5, 5), (:, :))       # both colons
    @test_throws ErrorException trim(ones(5, 5), (4:6, :))     # out of bounds
end

# ---------------------------------------------------------------------------
# crop / cropview
# ---------------------------------------------------------------------------
@testset "crop – correctness" begin
    x = reshape(1:25, 5, 5)

    @test crop(x, (3, 3)) == [7:9 12:14 17:19]
    @test crop(x, (1, 1)) == fill(13, 1, 1)
    @test crop(x, (4, 3), force_equal = false) == [6:9 11:14 16:19]

    # cropview returns a view
    @test cropview(x, (3, 3)) isa SubArray
    @test crop(x, (3, 3)) isa Array

    # Even-dimension frame
    y = reshape(1:16, 4, 4)
    @test crop(y, (2, 2)) == [6:7 10:11]

    # Float frame
    f = fill(2.0, 8, 8)
    @test crop(f, (4, 4)) == fill(2.0, 4, 4)

    # Warn on odd size adjustment
    @test_logs (:warn, "dimension 1 changed from 4 to 5") crop(x, (4, 3))

    # Errors
    @test_throws BoundsError       cropview(ones(5, 5), (7, 3))
    @test_throws ErrorException    cropview(ones(5, 5), (3, -5))
    @test_throws DimensionMismatch cropview(ones(5, 5), (3, 4, 5))
end

# ---------------------------------------------------------------------------
# combine
# ---------------------------------------------------------------------------
@testset "combine – correctness" begin
    frames = [reshape(1.0:4.0, 2, 2) for _ in 1:4]

    # median combine of identical frames = same frame
    @test combine(frames) ≈ [1.0 3.0; 2.0 4.0]

    # sum
    @test combine(frames, method = sum) ≈ [4.0 12.0; 8.0 16.0]

    # mean
    @test combine(frames, method = mean) ≈ [1.0 3.0; 2.0 4.0]

    # Vary frames and check median
    f1 = fill(1.0, 3, 3)
    f2 = fill(2.0, 3, 3)
    f3 = fill(3.0, 3, 3)
    @test combine([f1, f2, f3]) ≈ fill(2.0, 3, 3)

    # Vector vs varargs
    @test combine(f1, f2, f3) ≈ combine([f1, f2, f3])

    # Integer frames
    fi = [fill(Int32(k), 4, 4) for k in 1:5]
    @test combine(fi) ≈ fill(3.0, 4, 4)   # median of 1:5

    # Errors
    @test_throws DimensionMismatch combine(rand(5, 5), rand(6, 6))
    @test_throws MethodError combine()
end

# ---------------------------------------------------------------------------
# subtract_dark
# ---------------------------------------------------------------------------
@testset "subtract_dark – correctness" begin
    # Basic (same exposure)
    @test subtract_dark(ones(3, 3), ones(3, 3)) ≈ zeros(3, 3)

    # Scaled dark
    @test subtract_dark(ones(3, 3), ones(3, 3), data_exposure = 1, dark_exposure = 4) ≈ fill(0.75, 3, 3)
    @test subtract_dark(ones(Float32, 3, 3), ones(Float32, 3, 3), data_exposure = 13, dark_exposure = 17) ≈ fill(Float32(4/17), 3, 3)

    # Integer frame promoted to float
    fi = fill(Int32(10), 5, 5)
    di = fill(Int32(2), 5, 5)
    result = subtract_dark(fi, di, data_exposure = 2, dark_exposure = 1)
    @test result isa Matrix{Float64}
    @test result ≈ fill(6.0, 5, 5)

    # Mutating version
    f = fill(5.0, 3, 3)
    subtract_dark!(f, ones(3, 3), data_exposure = 2, dark_exposure = 1)
    @test f ≈ fill(3.0, 3, 3)

    # Type-stability
    @test @inferred(subtract_dark(ones(3, 3), ones(3, 3))) isa Matrix{Float64}
    @test @inferred(subtract_dark!(ones(3, 3), ones(3, 3))) isa Matrix{Float64}

    # Error: size mismatch
    @test_throws DimensionMismatch subtract_dark!(ones(5, 5), ones(6, 6))
end

# ---------------------------------------------------------------------------
# gain_correct
# ---------------------------------------------------------------------------
@testset "gain_correct – correctness" begin
    frame = fill(100.0, 4, 4)

    @test gain_correct(frame, 2.5) ≈ fill(250.0, 4, 4)
    @test gain_correct(frame, 1.0) ≈ frame

    # Integer frame -> float
    fi = fill(Int32(100), 4, 4)
    result = gain_correct(fi, 2.0)
    @test result isa Matrix{Float64}
    @test result ≈ fill(200.0, 4, 4)

    # Mutating version
    f = fill(10.0, 3, 3)
    gain_correct!(f, 3.0)
    @test f ≈ fill(30.0, 3, 3)

    # Non-mutating does not modify input
    orig = fill(5.0, 4, 4)
    _ = gain_correct(orig, 2.0)
    @test orig == fill(5.0, 4, 4)
end

# ---------------------------------------------------------------------------
# noise_model
# ---------------------------------------------------------------------------
@testset "noise_model – correctness" begin
    # Pure read noise (zero signal)
    sigma = noise_model(zeros(3, 3); read_noise = 10.0, gain = 1.0)
    @test sigma ≈ fill(10.0, 3, 3)

    # Pure shot noise (zero read noise, gain = 1)
    sigma2 = noise_model(fill(100.0, 3, 3); read_noise = 0.0, gain = 1.0)
    @test sigma2 ≈ fill(10.0, 3, 3)   # sqrt(100)

    # Combined
    v = noise_model(fill(100.0, 2, 2); read_noise = 10.0, gain = 2.0)
    # expected: sqrt(100/2 + (10/2)^2) = sqrt(50 + 25) = sqrt(75)
    @test v ≈ fill(sqrt(75.0), 2, 2)

    # Negative pixel values clamped to zero (no imaginary noise)
    sigma3 = noise_model(fill(-50.0, 2, 2); read_noise = 5.0, gain = 1.0)
    @test all(isfinite, sigma3)
    @test sigma3 ≈ fill(5.0, 2, 2)

    # Integer input
    sigma4 = noise_model(fill(Int32(400), 3, 3); read_noise = 0.0, gain = 1.0)
    @test sigma4 ≈ fill(20.0, 3, 3)

    # Output is always Float
    @test noise_model(ones(3, 3)) isa Matrix{Float64}
    @test noise_model(ones(Float32, 3, 3)) isa Matrix{Float32}
end

# ---------------------------------------------------------------------------
# cosmicray_lacosmic
# ---------------------------------------------------------------------------
@testset "cosmicray_lacosmic – correctness" begin
    # Single obvious cosmic ray
    frame = fill(100.0, 30, 30)
    frame[15, 15] = 50_000.0

    cleaned, mask = cosmicray_lacosmic(frame; gain = 1.0, read_noise = 5.0, sigma_clip = 4.5)

    @test mask[15, 15]                 # cosmic ray detected
    @test cleaned[15, 15] < 50_000.0  # pixel replaced
    # Background pixels should not be flagged
    @test sum(mask) < 10

    # Multiple cosmic rays
    frame2 = fill(200.0, 40, 40)
    frame2[10, 10] = 1e5
    frame2[20, 20] = 1e5
    frame2[30, 30] = 1e5

    cleaned2, mask2 = cosmicray_lacosmic(frame2; gain = 2.0, read_noise = 8.0)
    @test mask2[10, 10]
    @test mask2[20, 20]
    @test mask2[30, 30]
    @test all(cleaned2 .< 1e5)

    # Clean frame: no cosmic rays should be flagged
    flat_frame = fill(500.0, 20, 20)
    _, mask3 = cosmicray_lacosmic(flat_frame; gain = 1.0, read_noise = 5.0)
    @test !any(mask3)

    # Output type matches float(eltype(input))
    frame_f32 = fill(Float32(100), 20, 20)
    frame_f32[10, 10] = Float32(1e5)
    cleaned_f32, _ = cosmicray_lacosmic(frame_f32)
    @test cleaned_f32 isa Matrix{Float32}
end

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
@testset "helper functions" begin
    # axes_min_length
    @test axes_min_length((:, :)) == 1
    @test axes_min_length((:, 5:10)) == 2
    @test axes_min_length((1:2, :)) == 1
    @test axes_min_length((1:2, 5:12)) == 1
    @test axes_min_length((10:100, 1:5)) == 2

    # fits_indices
    @test fits_indices("[1024:2048, 200:300]") == [200:300, 1024:2048]
    @test fits_indices("[:, 200:300]") == [200:300, :]
    @test fits_indices("[1024:2048, :]") == [:, 1024:2048]
    @test fits_indices("[200:300, 1024:2048]") == [1024:2048, 200:300]

    # convert_value
    @test convert_value(Int16, 5.4) == 5
    @test convert_value(Int32, 5.4) == 5
    @test convert_value(Int64, -5.4) == -5
    @test convert_value(Float64, -5.4) ≈ -5.4
    @test convert_value(Float32, -5.4) ≈ -5.4
end
