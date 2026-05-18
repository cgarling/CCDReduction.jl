using FITSIO

# ---------------------------------------------------------------------------
# Generate mock FITS files used across the test suite
# ---------------------------------------------------------------------------
# Instead of downloading real data from the internet we write small synthetic
# FITS files to a temporary directory that persists for the lifetime of the
# test run.

const TEST_DATA_DIR = mktempdir()

"""
    make_fits(path, data; keywords...)

Write a 2-D `data` array to a new single-HDU FITS file at `path`, inserting
any supplied keyword-value pairs into the primary header.
"""
function make_fits(path::String, data::AbstractMatrix; kwargs...)
    header = FITSHeader(
        ["SIMPLE", string.(keys(kwargs))...],
        [true,     values(kwargs)...],
        ["",       fill("", length(kwargs))...],
    )
    FITS(path, "w") do fh
        d = ndims(data)
        write(fh, permutedims(data, d:-1:1); header = header)
    end
    return path
end

# --- Frame A: 64×64 Float32 image with known constant value ----------------
const FRAME_A_DATA = fill(Float32(500.0), 64, 64)
const FRAME_A_PATH = make_fits(
    joinpath(TEST_DATA_DIR, "frameA.fits"), FRAME_A_DATA;
    EXPTIME = 60.0, IMAGETYP = "LIGHT", FILTER = "V",
)

# --- Frame B: 64×64 Int32 image (photon-count style) -----------------------
const FRAME_B_DATA = fill(Int32(1000), 64, 64)
const FRAME_B_PATH = make_fits(
    joinpath(TEST_DATA_DIR, "frameB.fits"), FRAME_B_DATA;
    EXPTIME = 120.0, IMAGETYP = "LIGHT", FILTER = "R",
)

# --- Bias frame: 64×64 zeros -----------------------------------------------
const BIAS_DATA = zeros(Float32, 64, 64)
const BIAS_PATH = make_fits(
    joinpath(TEST_DATA_DIR, "bias.fits"), BIAS_DATA;
    EXPTIME = 0.0, IMAGETYP = "BIAS",
)

# --- Dark frame: uniform 2.0 e⁻/s ------------------------------------------
const DARK_DATA = fill(Float32(2.0), 64, 64)
const DARK_PATH = make_fits(
    joinpath(TEST_DATA_DIR, "dark.fits"), DARK_DATA;
    EXPTIME = 1.0, IMAGETYP = "DARK",
)

# --- Flat frame: uniform 1.0 (already normalised) --------------------------
const FLAT_DATA = fill(Float32(1.0), 64, 64)
const FLAT_PATH = make_fits(
    joinpath(TEST_DATA_DIR, "flat.fits"), FLAT_DATA;
    EXPTIME = 30.0, IMAGETYP = "FLAT",
)
