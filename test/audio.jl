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
    audio_one_channel = SampleBuf(rand(10000), 1000)
    audio_two_channel = SampleBuf(rand(10000,2), 1000)
    audio_multi_channel = SampleBuf(rand(10000,2), 1000)
    # Tests for pitchshift method
    @testset "pitchshift" begin
        @test isapprox(duration(pitchshift(audio_one_channel)), duration(audio_one_channel), atol=0.2)
        @test isapprox(duration(pitchshift(audio_two_channel)), duration(audio_two_channel), atol=0.2)
        @test isapprox(duration(pitchshift(audio_multi_channel)), duration(audio_multi_channel), atol=0.2)
    end

    audio_one_channel = SampleBuf(rand(10000), 1000)
    audio_two_channel = SampleBuf(rand(10000,2), 1000)
    audio_multi_channel = SampleBuf(rand(10000,2), 1000)
    # Tests for speedup method
    @testset "speedup" begin
        
        @test isapprox(duration(speedup(audio_one_channel, 2)), duration(audio_one_channel) / 2.0 , atol=0.4)
        @test isapprox(duration(speedup(audio_two_channel, 2)), duration(audio_two_channel) / 2.0 , atol=0.4)
        @test isapprox(duration(speedup(audio_multi_channel, 2)), duration(audio_multi_channel) / 2.0 , atol=0.4)
    end

    audio_one_channel = SampleBuf(rand(10000), 1000)
    audio_two_channel = SampleBuf(rand(10000,2), 1000)
    audio_multi_channel = SampleBuf(rand(10000,2), 1000)
    # Tests for slowdown method
    @testset "slowdown" begin

        @test isapprox(duration(slowdown(audio_one_channel, 2)), duration(audio_one_channel) * 2.0, atol=1.2)
        @test isapprox(duration(slowdown(audio_two_channel, 2)), duration(audio_two_channel) * 2.0, atol=1.2)
        @test isapprox(duration(slowdown(audio_multi_channel, 2)), duration(audio_multi_channel) * 2.0, atol=1.2)
   end

    # Tests for zero crossing rate method
    @testset "zero_crossing_rate" begin
        audio_one_channel = SampleBuf(rand(100000), 10)

        expected = length(zero_crossing_rate(audio_one_channel, 2048, 512))
        res = 192
        @test expected == res
    end

end