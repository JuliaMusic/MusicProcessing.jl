module MusicProcessing

using DSP
using FFTW
using FixedPointNumbers
using Requires
using SampledSignals

# types used for fixed-point 16-bit and 32-bit encoding
const PCM16Sample = Fixed{Int16, 15}
const PCM32Sample = Fixed{Int32, 31}

export Hz, kHz, s, ..
export PCM16Sample, PCM32Sample

include("util.jl")
include("complex.jl")

include("audio.jl")
include("TFR.jl")

include("mel.jl")
include("constantq.jl")
include("chroma.jl")

function __init__()
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("display.jl")
end

end # module
