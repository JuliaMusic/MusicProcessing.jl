# audio.jl

@testset verbose= true "audio.jl" begin
    
    # Example audio
    audio_one_channel = SampleBuf(rand(1000), 10)
    audio_two_channel = SampleBuf(rand(1000,2), 10)
    audio_multi_channel = SampleBuf(rand(1000,2), 10)
    
    # Tests for mono
    @testset "mono" begin
        @test size(mono(audio_two_channel)) == (1000,)
        @test size(mono(audio_multi_channel)) == (1000,)
        
    end

    # Tests for resample method
    @testset "resample.jl" begin
        @test size(resample(audio_one_channel, 5))[1] == 500
        @test size(resample(audio_two_channel, 5))[1] == 500
        @test size(resample(audio_multi_channel, 5))[1] == 500
    end

    # Tests for duration method
    @testset "duration" begin
        @test duration(audio_one_channel) == 100
        @test duration(audio_two_channel) == 100
        @test duration(audio_multi_channel) == 100
    end

    # Tests for pitchshift method
    @testset "pitchshift" begin
        @test duration(pitchshift(audio_one_channel)) == duration(audio_one_channel)
        @test duration(pitchshift(audio_two_channel)) == duration(audio_two_channel)
        @test duration(pitchshift(audio_multi_channel)) == duration(audio_multi_channel)
    end

    # Tests for speedup method
    @testset "speedup" begin
        @test duration(speedup(audio_one_channel, 2)) ≈ duration(audio_one_channel) // 2
        @test duration(speedup(audio_two_channel, 2)) ≈ duration(audio_two_channel) // 2
        @test duration(speedup(audio_multi_channel, 2)) ≈ duration(audio_multi_channel) // 2
    end

    # Tests for slowdown method
    @testset "slowdown" begin
        @test duration(slowdown(audio_one_channel, 2)) ≈ duration(audio_one_channel) * 2
        @test duration(slowdown(audio_two_channel, 2)) ≈ duration(audio_two_channel) * 2
        @test duration(slowdown(audio_multi_channel, 2)) ≈ duration(audio_multi_channel) * 2
    end

    # Tests for zero crossing rate method
    @testset "zero_crossing_rate" begin
        audio_one_channel = SampleBuf(rand(100000), 10)

        expected = length(zero_crossing_rate(audio_one_channel, 2048, 512))
        res = 192
        @test expected == res
    end
end