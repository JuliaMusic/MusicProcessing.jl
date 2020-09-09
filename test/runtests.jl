using MusicProcessing
using SampledSignals: SampleBuf
using Test

# Example audio
audio_one_channel = SampleBuf(rand(1000), 10)
audio_two_channel = SampleBuf(rand(1000,2), 10)
audio_multi_channel = SampleBuf(rand(1000,2), 10)

# audio.jl

@testset "audio.jl" begin
    @test size(mono(audio_two_channel)) == (1000,)
    @test size(mono(audio_multi_channel)) == (1000,)

    @test size(resample(audio_one_channel, 5))[1] == 500
    @test size(resample(audio_two_channel, 5))[1] == 500
    @test size(resample(audio_multi_channel, 5))[1] == 500

    @test duration(audio_one_channel) == 100
    @test duration(audio_two_channel) == 100
    @test duration(audio_multi_channel) == 100

    @test duration(pitchshift(audio_one_channel)) == duration(audio_one_channel)
    @test duration(pitchshift(audio_two_channel)) == duration(audio_two_channel)
    @test duration(pitchshift(audio_multi_channel)) == duration(audio_multi_channel)

    @test duration(speedup(audio_one_channel, 2)) ≈ duration(audio_one_channel) // 2
    @test duration(speedup(audio_two_channel, 2)) ≈ duration(audio_two_channel) // 2
    @test duration(speedup(audio_multi_channel, 2)) ≈ duration(audio_multi_channel) // 2

    @test duration(slowdown(audio_one_channel, 2)) ≈ duration(audio_one_channel) * 2
    @test duration(slowdown(audio_two_channel, 2)) ≈ duration(audio_two_channel) * 2
    @test duration(slowdown(audio_multi_channel, 2)) ≈ duration(audio_multi_channel) * 2
end
