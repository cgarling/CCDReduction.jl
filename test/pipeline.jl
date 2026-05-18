# Tests for Pipeline and ChunkedPipeline

@testset "Pipeline – basics" begin
    bias = fill(1.0, 8, 8)
    flat = fill(2.0, 8, 8)

    p = Pipeline(
        frame -> subtract_bias(frame, bias),
        frame -> flat_correct(frame, flat; norm_value = 1.0),
    )

    @test length(p) == 2

    science = fill(3.0, 8, 8)
    # bias-subtracted → 2.0; flat-correct with norm=1,flat=2 → 2/2=1
    result = process(p, science)
    @test result ≈ fill(1.0, 8, 8)

    # Callable syntax
    result2 = p(science)
    @test result2 ≈ result
end

@testset "Pipeline – identity" begin
    p = Pipeline(identity)
    x = rand(5, 5)
    @test process(p, x) === x
end

@testset "Pipeline – ordering" begin
    # Step 1: add 10; step 2: multiply by 2 → (x+10)*2 ≠ (x*2)+10 in general
    p = Pipeline(x -> x .+ 10, x -> x .* 2)
    frame = fill(1.0, 3, 3)
    @test process(p, frame) ≈ fill(22.0, 3, 3)  # (1+10)*2 = 22
end

@testset "ChunkedPipeline – basics" begin
    bias = zeros(64, 64)
    flat = fill(2.0, 64, 64)

    p = ChunkedPipeline(
        tiled_step(subtract_bias, bias),
        tiled_step(flat_correct, flat);
        chunksize = (16, 16),
    )

    @test length(p) == 2

    science = fill(4.0, 64, 64)
    result = process(p, science)
    # bias=0 → 4; flat_correct with mean(flat)=2 → 4/(2/2)=4
    @test result ≈ fill(4.0, 64, 64)
end

@testset "ChunkedPipeline – matches Pipeline output" begin
    science = rand(64, 64)
    bias    = rand(64, 64) .* 0.1
    # Use a constant flat so that mean(flat[ri,rj]) == mean(flat) for any region
    flat    = fill(3.0, 64, 64)

    sp = Pipeline(
        frame -> subtract_bias(frame, bias),
        frame -> flat_correct(frame, flat),
    )

    cp = ChunkedPipeline(
        tiled_step(subtract_bias, bias),
        tiled_step(flat_correct, flat);
        chunksize = (16, 16),
    )

    result_seq   = process(sp, science)
    result_chunk = process(cp, science)
    @test result_chunk ≈ result_seq  atol=1e-10
end

@testset "ChunkedPipeline – non-square and non-divisible size" begin
    science = rand(30, 50)
    bias    = zeros(30, 50)

    cp = ChunkedPipeline(
        tiled_step(subtract_bias, bias);
        chunksize = (7, 11),  # deliberately non-divisible
    )

    result = process(cp, science)
    @test size(result) == (30, 50)
    @test result ≈ science
end

@testset "ChunkedPipeline – output is float" begin
    # Int32 input: float(Int32) == Float64
    science_int = fill(Int32(100), 16, 16)
    p = ChunkedPipeline(chunk -> chunk; chunksize = (4, 4))
    result = process(p, science_int)
    @test result isa Matrix{Float64}
    @test all(result .== 100.0)
end

@testset "tiled_step" begin
    bias = zeros(32, 32)
    flat = fill(2.0, 32, 32)
    science = fill(4.0, 32, 32)

    p = ChunkedPipeline(
        tiled_step(subtract_bias, bias),
        tiled_step(flat_correct, flat);
        chunksize = (8, 8),
    )

    result = process(p, science)
    # subtract_bias: 4-0=4; flat_correct with mean(flat)=2 → 4/(2/2)=4
    @test result ≈ fill(4.0, 32, 32)
end

@testset "Pipeline – CCDData passthrough" begin
    bias = CCDData(zeros(8, 8))
    p = Pipeline(frame -> subtract_bias(frame, bias))
    science = CCDData(fill(5.0, 8, 8))
    result = process(p, science)
    @test result isa CCDData
    @test result.data ≈ fill(5.0, 8, 8)
end

@testset "ChunkedPipeline – error on non-positive chunksize" begin
    @test_throws ArgumentError ChunkedPipeline(identity; chunksize = (0, 16))
    @test_throws ArgumentError ChunkedPipeline(identity; chunksize = (16, -1))
end
