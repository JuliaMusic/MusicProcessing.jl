# simple methods to process audio signals

import DSP.resample, DSP.arraysplit
export mono, resample, duration, play, pitchshift, speedup, slowdown, zero_crossing_rate

"""
    mono(audio)

convert a multichannel audio to mono
"""
function mono(audio::SampleBuf{T, 2}) where {T <: AbstractFloat}
    SampleBuf{T, 1}(
        vec(monoch(audio).data),
        audio.samplerate
    )
end

# special code for fixed-point samples, to avoid overflow
function mono(audio::SampleBuf{T, 2}) where {T <: Fixed}
    nchannels = nchannels(audio)
    if nchannels == 1
        SampleBuf{T, 1}(vec(audio.data), audio.samplerate)
    elseif nchannels == 2
        nsamples = nframes(audio)
        buffer = Array{T}(undef, nsamples)
        for i = 1:nsamples
            @inbounds a = audio.data[i, 1].i
            @inbounds b = audio.data[i, 2].i
            m = (a >> 1) + (b >> 1) + (a & b & 1)
            @inbounds buffer[i] = T(m, 0)
        end
        SampleBuf{T, 1}(buffer, audio.samplerate)
    else
        SampleBuf{T, 1}(
            map(T, vec(mean(map(Float32, audio.data)))),
            audio.samplerate
        )
    end
end

"""resample audio with a different sample rate"""
function resample(audio::SampleBuf{T, 2}, samplerate::Real) where T
    rate = samplerate / audio.samplerate
    filter = DSP.resample_filter(rate, 512, 1.0, 140)
    SampleBuf{T, 2}(
        mapslices(audio.data, dims=1) do data
            DSP.resample(data, rate, filter)
        end,
        samplerate
    )
end

"""resample audio with a different sample rate"""
function resample(audio::SampleBuf{T, 1}, samplerate::Real) where {T, F}
    SampleBuf{T, 1}(
        DSP.resample(audio.data, samplerate / audio.samplerate),
        samplerate
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
    nchannels = nchannels(audio)
    stream = PortAudio.PortAudioStream(2, nchannels)
    try
        write(stream, audio)
    finally
        close(stream)
    end
end
play(audio::SampleBuf{T}) where T = play(map(Float32, audio))


""""""
function pitchshift(audio::SampleBuf{T, N}, semitones::Real = 8) where {T, N}
    rate = 2.0 ^ (semitones / 12.0)
    shifted = resample(slowdown(audio, rate), audio.samplerate / rate)
    SampleBuf{T, N}(
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
function zero_crossing_rate(audio::SampleBuf{T, 1}, framesize::Int = 1024, hopsize::Int = framesize >> 2) where T
    nframes = MusicProcessing.nframes(length(audio.data), framesize, hopsize)
    result = Array{Float32}(undef,nframes)

    offset = 0
    for i = 1:nframes
        result[i] = zero_crossings(audio.data, framesize, offset) / framesize
        offset += hopsize
    end

    result
end
