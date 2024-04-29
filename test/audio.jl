# audio.jl
using MusicProcessing, Test

# Example audio
audio_one   = SampleBuf(rand(1000  ), 10)
audio_two   = SampleBuf(rand(1000,2), 10)
audio_multi = SampleBuf(rand(1000,3), 10)

@testset "mono" begin
    for audio in [audio_one, audio_two, audio_multi]
        @test size(mono(audio)) == (1000,)
    end
end

@testset "resample" begin
    for audio in [audio_one, audio_two, audio_multi]
        @test size(resample(audio, 5), 1) == 500
    end
end

@testset "duration" begin
    for audio in [audio_one, audio_two, audio_multi]
        @test duration(audio) == 100
    end
end

# Example audio
audio_one   = SampleBuf(rand(10000  ), 1000)
audio_two   = SampleBuf(rand(10000,2), 1000)
audio_multi = SampleBuf(rand(10000,3), 1000)

@testset "pitchshift" begin
    for audio in [audio_one, audio_two, audio_multi]
        @test isapprox(duration(pitchshift(audio)), duration(audio), atol=0.2)
    end
end

@testset "speedup" begin
    for audio in [audio_one, audio_two, audio_multi]
        @test isapprox(duration(speedup(audio, 2)), duration(audio) / 2.0 , atol=0.4)
    end
end

@testset "slowdown" begin
    for audio in [audio_one, audio_two, audio_multi]
        @test isapprox(duration(slowdown(audio, 2)), duration(audio) * 2.0, atol=1.2)
    end
end

# Example audio
audio_one = SampleBuf(rand(100000), 10)

@testset "zero_crossing_rate" begin
    @test length(zero_crossing_rate(audio_one, 2048, 512)) == 192
end
