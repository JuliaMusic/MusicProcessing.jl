module MusicProcessingTest

using MusicProcessing
using Test
using WAV

tests = [
    "audio.jl",
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end

end


