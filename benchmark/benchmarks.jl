using BenchmarkTools
using CCDReduction
using Statistics

# AirspeedVelocity.jl-compatible benchmark suite.
# Top-level keys correspond to function / type names; sub-keys encode
# the frame size and element type so that individual benchmarks are
# uniquely identified.

const SUITE = BenchmarkGroup()

# Pre-create top-level groups
for grp in ("subtract_bias", "subtract_dark", "flat_correct", "gain_correct",
            "noise_model", "subtract_overscan", "trim", "combine",
            "Pipeline", "ChunkedPipeline", "cosmicray_lacosmic")
    SUITE[grp] = BenchmarkGroup()
end

# ---------------------------------------------------------------------------
# Frame sizes
# ---------------------------------------------------------------------------
const SIZES = (("512x512", (512, 512)), ("2048x2048", (2048, 2048)))

for (label, sz) in SIZES
    nrows, ncols = sz

    frame_f  = ones(Float64, sz...)
    frame_i  = ones(Int32,   sz...)
    bias_f   = fill(100.0,   sz...)
    dark_f   = fill(  2.0,   sz...)
    flat_f   = fill(  3.0,   sz...)
    gain_val = 2.0

    # --- subtract_bias ------------------------------------------------------
    sg = SUITE["subtract_bias"]
    sg["float_$label"]   = @benchmarkable subtract_bias($frame_f, $bias_f)
    sg["int_$label"]     = @benchmarkable subtract_bias($frame_i, $bias_f)
    sg["inplace_$label"] = @benchmarkable subtract_bias!(x, $bias_f) setup=(x=copy($frame_f))

    # --- subtract_dark ------------------------------------------------------
    sg = SUITE["subtract_dark"]
    sg["float_$label"]  = @benchmarkable subtract_dark($frame_f, $dark_f)
    sg["scaled_$label"] = @benchmarkable subtract_dark($frame_f, $dark_f,
                                                        data_exposure=60.0,
                                                        dark_exposure=300.0)

    # --- flat_correct -------------------------------------------------------
    sg = SUITE["flat_correct"]
    sg["float_$label"] = @benchmarkable flat_correct($frame_f, $flat_f)
    sg["int_$label"]   = @benchmarkable flat_correct($frame_i, $flat_f)

    # --- gain_correct -------------------------------------------------------
    sg = SUITE["gain_correct"]
    sg["float_$label"] = @benchmarkable gain_correct($frame_f, $gain_val)
    sg["int_$label"]   = @benchmarkable gain_correct($frame_i, $gain_val)

    # --- noise_model --------------------------------------------------------
    sg = SUITE["noise_model"]
    sg["float_$label"] = @benchmarkable noise_model($frame_f; read_noise=6.5, gain=1.5)

    # --- subtract_overscan --------------------------------------------------
    frame_ov = ones(Float64, nrows, ncols + 20)
    sg = SUITE["subtract_overscan"]
    sg["float_$label"] = @benchmarkable subtract_overscan($frame_ov,
                                                           (:, $(ncols+1:ncols+20)))

    # --- trim ---------------------------------------------------------------
    frame_tr = ones(Float64, nrows, ncols + 20)
    sg = SUITE["trim"]
    sg["float_$label"] = @benchmarkable trim($frame_tr, (:, $(ncols+1:ncols+20)))

    # --- combine (3 frames, median / mean / sum) ----------------------------
    frames3 = [copy(frame_f) for _ in 1:3]
    sg = SUITE["combine"]
    sg["median_3x_$label"] = @benchmarkable combine($frames3)
    sg["mean_3x_$label"]   = @benchmarkable combine($frames3, method=mean)
    sg["sum_3x_$label"]    = @benchmarkable combine($frames3, method=sum)

    # --- Pipeline (bias + dark + flat) --------------------------------------
    pipeline = Pipeline(
        f -> subtract_bias(f, bias_f),
        f -> subtract_dark(f, dark_f; data_exposure=60.0, dark_exposure=300.0),
        f -> flat_correct(f, flat_f),
    )
    SUITE["Pipeline"]["full_reduction_$label"] = @benchmarkable process($pipeline, $frame_f)

    # --- ChunkedPipeline (same steps, tiled 256×256) ------------------------
    chunked = ChunkedPipeline(
        tiled_step(subtract_bias, bias_f),
        tiled_step(subtract_dark, dark_f),
        tiled_step(flat_correct,  flat_f);
        chunksize = (256, 256),
    )
    SUITE["ChunkedPipeline"]["full_reduction_$label"] = @benchmarkable process($chunked, $frame_f)
end

# ---------------------------------------------------------------------------
# cosmicray_lacosmic (expensive; only 256×256)
# ---------------------------------------------------------------------------
let sz = (256, 256)
    frame_cr = fill(500.0, sz...)
    frame_cr[128, 128] = 1e6
    SUITE["cosmicray_lacosmic"]["256x256"] =
        @benchmarkable cosmicray_lacosmic($frame_cr; gain=1.5, read_noise=6.5, niter=4)
end
