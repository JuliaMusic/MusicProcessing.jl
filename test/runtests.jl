module MusicProcessingTest

using MusicProcessing
using SampledSignals: SampleBuf
using Test

tests = [
    "audio.jl",
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end

end


