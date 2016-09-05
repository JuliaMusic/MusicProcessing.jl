# simple methods to process audio signals

import DSP.resample, DSP.arraysplit
export mono, resample, duration, play, pitchshift, speedup, slowdown, zero_crossing_rate

"""
    mono(audio)

convert a multichannel audio to mono
"""
function mono{T <: AbstractFloat}(audio::SampleBuf{T, 2, Hertz})
    SampleBuf{T, 1, Hertz}(
        mean(audio.data, 2)[:],
        audio.samplerate
    )
end

# special code for fixed-point samples, to avoid overflow
function mono{T <: Fixed}(audio::SampleBuf{T, 2, Hertz})
    nchannels = SampledSignals.nchannels(audio)
    if nchannels == 1
        SampleBuf{T, 1, Hertz}(audio.data[:], audio.samplerate)
    elseif nchannels == 2
        nsamples = SampledSignals.nframes(audio)
        buffer = Array(T, nsamples)
        for i = 1:nsamples
            @inbounds a = audio.data[i, 1].i
            @inbounds b = audio.data[i, 2].i
            m = (a >> 1) + (b >> 1) + (a & b & 1)
            @inbounds buffer[i] = T(m, 0)
        end
        SampleBuf{T, 1, Hertz}(buffer, audio.samplerate)
    else
        SampleBuf{T, 1, Hertz}(
            map(T, mean(map(Float32, audio.data))[:]),
            audio.samplerate
        )
    end
end

"""resample audio with a different sample rate"""
function resample{T, F}(audio::SampleBuf{T, 2, Hertz}, samplerate::SIUnits.SIQuantity{F,0,0,-1,0,0,0,0,0,0})
    sr = hertz(float(samplerate))
    rate = sr / audio.samplerate
    filter = DSP.resample_filter(rate, 512, 1.0, 140)
    SampleBuf{T, 2, SIUnits.SIQuantity{F,0,0,-1,0,0,0,0,0,0}}(
        mapslices(audio.data, 1) do data
            DSP.resample(data, rate, filter)
        end,
        sr
    )
end

"""resample audio with a different sample rate"""
function resample{T, F}(audio::SampleBuf{T, 1, Hertz}, samplerate::SIUnits.SIQuantity{F,0,0,-1,0,0,0,0,0,0})
    sr = hertz(float(samplerate))
    SampleBuf{T, 1, SIUnits.SIQuantity{F,0,0,-1,0,0,0,0,0,0}}(
        DSP.resample(audio.data, sr / audio.samplerate),
        sr
    )
end

"""returns the duration of given audio, in seconds"""
function duration(audio::SampleBuf)
    nframes(audio) / samplerate(audio)
end

"""
    play(audio)

play the audio on local computer using PortAudio
"""
function play(audio::SampleBuf{Float32})
    # import PortAudio on-demand
    @eval import PortAudio
    nchannels = SampledSignals.nchannels(audio)
    stream = PortAudio.PortAudioStream(2, nchannels)
    try
        write(stream, audio)
    finally
        close(stream)
    end
end
play{T}(audio::SampleBuf{T}) = play(map(Float32, audio))


""""""
function pitchshift{T, N}(audio::SampleBuf{T, N, Hertz}, semitones::Real)
    rate = 2.0 ^ (semitones / 12.0)
    shifted = resample(slowdown(audio, rate), audio.samplerate / rate)
    SampleBuf{T, N, Hertz}(
        shifted.data,
        audio.samplerate
    )
end

""""""
function speedup(audio::SampleBuf, speed::Real, windowsize::Int = 1024, hopsize::Int = windowsize >> 2; kwargs...)
    S = stft(audio, windowsize, hopsize; kwargs...)
    S = phase_vocoder(S, speed, hopsize)
    istft(S, audio.samplerate, windowsize, hopsize; kwargs...)
end

""""""
function slowdown(audio::SampleBuf, ratio::Real, windowsize::Int = 1024, hopsize::Int = windowsize >> 2; kwargs...)
    speedup(audio, 1.0 / ratio, windowsize, hopsize; kwargs...)
end

""""""
function zero_crossing_rate{T}(audio::SampleBuf{T, 1}, framesize::Int = 1024, hopsize::Int = framesize >> 2)
    nframes = MusicProcessing.nframes(length(audio.data), framesize, hopsize)
    result = Array(Float32, nframes)

    offset = 0
    for i = 1:nframes
        result[i] = zero_crossings(audio.data, framesize, offset) / framesize
        offset += hopsize
    end

    SampleBuf{T, 1, SIUnits.SIQuantity{Float32,0,0,-1}}(result, audio.samplerate / Float32(hopsize))
end
