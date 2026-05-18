module CCDReduction

using Statistics
using LazyStack
using FITSIO

export subtract_bias,
       subtract_bias!,
       subtract_overscan,
       subtract_overscan!,
       flat_correct,
       flat_correct!,
       trim,
       trimview,
       crop,
       cropview,
       combine,
       subtract_dark,
       subtract_dark!,
       gain_correct,
       gain_correct!,
       noise_model,
       cosmicray_lacosmic,
       fitscollection,
       arrays,
       filenames,
       ccds,
       CCDData,
       data,
       hdr,
       ImageCollection,
       Pipeline,
       ChunkedPipeline,
       tiled_step,
       process

include("ccddata.jl")
include("methods.jl")
include("fits.jl")
include("collection.jl")
include("pipeline.jl")

end
