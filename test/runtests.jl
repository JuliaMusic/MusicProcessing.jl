using MusicProcessing
using Base.Test

# write your own tests here
@test 1 == 1



include("MusicProcessing.jl")
using MusicProcessing
using SampledSignals

# Example Vars

audio = SampleBuf(rand(10), 0.5)
# TFR
spectrogram(audio)
