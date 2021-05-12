# simple methods to process audio signals

import DSP.resample, DSP.arraysplit


"""
    mono(audio::SampleBuf{T, N})

Converts a multichannel audio to single channel

# Output 

Return the mono channel audio as a SampleBuf{T, 1} with initial sample rate.

# Arguments

## `audio`

An audio of type SampleBuf{T, N} which is to be converted to single channel audio.

# Example

```julia
using MusicProcessing

audio_one_channel = SampleBuf(rand(1000), 10)
audio_two_channel = SampleBuf(rand(1000, 2), 10)
audio_multi_channel = SampleBuf(rand(1000, 4), 10)

mono(audio_one_channel)
mono(audio_two_channel)
mono(audio_multi_channel)
```

"""
function mono(audio::SampleBuf{T, N}) where {T <: AbstractFloat, N}
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

"""
    resample(audio::SampleBuf{T, 2}, samplerate::Real) 

Resamples input audio with a provided sample rate

# Output 

Return the resampled audio as a SampleBuf{T, 2} with given sample rate.

# Arguments

## `audio`

An audio of type SampleBuf{T, 2} which is to be resampled.

## `samplerate`

A parameter to specify the samplerate for final output's sample rate

# Example

```julia
using MusicProcessing

audio_two_channel = SampleBuf(rand(1000, 2), 10)

resample(audio_two_channel, 5)
```
"""
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

"""
    resample(audio::SampleBuf{T, 2}, samplerate::Real) 

Resamples input audio with a provided sample rate

"""
function resample(audio::SampleBuf{T, 1}, samplerate::Real) where {T, F}
    SampleBuf{T, 1}(
        DSP.resample(audio.data, samplerate / audio.samplerate),
        samplerate
    )
end

"""
    duration(audio::SampleBuf)

Returns the duration of given audio, in seconds

# Output 

Return the duration of audio(in seconds).

# Arguments

## `audio`

An audio of type SampleBuf whose duration is to be found.

# Example

```julia
using MusicProcessing

audio_one_channel = SampleBuf(rand(1000), 10)
audio_two_channel = SampleBuf(rand(1000, 2), 10)
audio_multi_channel = SampleBuf(rand(1000, 4), 10)

duration(audio_one_channel)
duration(audio_two_channel)
duration(audio_multi_channel)
```

"""
function duration(audio::SampleBuf)
    nframes(audio) / samplerate(audio)
end

"""
    play(audio::SampleBuf{Float32})

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


"""
    pitchshift(audio::SampleBuf{T, N}, semitones::Real = 8)

Returns pitchshifted audio waveform

# Output 

Return the pitchshifted audio waveform as a SampleBuf{T, N}.

# Arguments

## `audio`

An audio of type SampleBuf{T, N} whose pitch is to be shifted.

## `semitones`

A parameter used to define the new sample rate for the audio.\n
`new samplerate` = 2.0 ^ (semitones / 12.0)

# Example

```julia
using MusicProcessing

audio_one_channel = SampleBuf(rand(10000), 10)
audio_two_channel = SampleBuf(rand(10000, 2), 10)
audio_multi_channel = SampleBuf(rand(10000, 4), 10)

pitchshift(audio_one_channel, 8)
pitchshift(audio_two_channel, 8)
pitchshift(audio_multi_channel, 8)

```

"""
function pitchshift(audio::SampleBuf{T, N}, semitones::Real = 8) where {T, N}
    rate = 2.0 ^ (semitones / 12.0)
    shifted = resample(slowdown(audio, rate), audio.samplerate / rate)
    SampleBuf{T, N}(
        shifted.data,
        audio.samplerate
    )
end

"""
    speedup(audio::SampleBuf, speed::Real, windowsize::Int = 1024, hopsize::Int = windowsize >> 2; kwargs...)

Returns a speedup audio waveform using various parameters

# Output 

Return the speedup audio waveform as a SampleBuf{T, N} using windowsize and hopsize.

# Arguments

## `audio`

An audio of type SampleBuf{T, N} whose pitch is to be shifted.

## `windowsize`

A parameter used to define the window size for stft and istft functions.


## `hopsize`

A parameter used to define the number audio of frames between stft columns. 

# Example

```julia
using MusicProcessing

audio_one_channel = SampleBuf(rand(10000), 10)
audio_two_channel = SampleBuf(rand(10000, 2), 10)
audio_multi_channel = SampleBuf(rand(10000, 4), 10)

speedup(audio_one_channel, 2048, 512)
speedup(audio_one_channel, 2048, 512)
speedup(audio_one_channel, 2048, 512)

```

"""
function speedup(audio::SampleBuf, speed::Real, windowsize::Int = 1024, hopsize::Int = windowsize >> 2; kwargs...)
    S = stft(audio, windowsize, hopsize; kwargs...)
    S = phase_vocoder(S, speed, hopsize)
    istft(S, audio.samplerate, windowsize, hopsize; kwargs...)
end

"""
    slowdown(audio::SampleBuf, ratio::Real, windowsize::Int = 1024, hopsize::Int = windowsize >> 2; kwargs...)

Returns a slowdowned audio using various parameters

# Output 

Return the slowdowned audio waveform as a SampleBuf{T, N} using windowsize and hopsize.

# Arguments

## `audio`

An audio of type SampleBuf{T, N} whose speed is to be slowdowned.

## `windowsize`

A parameter used to define the window size for stft and istft functions.


## `hopsize`

A parameter used to define the number audio of frames between stft columns. 

# Example

```julia
using MusicProcessing

audio_one_channel = SampleBuf(rand(10000), 10)
audio_two_channel = SampleBuf(rand(10000, 2), 10)
audio_multi_channel = SampleBuf(rand(10000, 4), 10)

slowdown(audio_one_channel, 2048, 512)
slowdown(audio_two_channel, 2048, 512)
slowdown(audio_multi_channel, 2048, 512)

```

"""
function slowdown(audio::SampleBuf, ratio::Real, windowsize::Int = 1024, hopsize::Int = windowsize >> 2; kwargs...)
    speedup(audio, 1.0 / ratio, windowsize, hopsize; kwargs...)
end

"""
    zero_crossing_rate(audio::SampleBuf{T, 1}, framesize::Int = 1024, hopsize::Int = framesize >> 2) 

Returns zero crossing rate in a audio waveform

# Output 

Returns a array of zero-crossings in y along the selected axis of audio using various parameters

# Arguments

## `audio`

An audio of type SampleBuf{T, N} whose speed is to be slowdowned.

## `framesize`

A parameter used to define the framesize to calculate no of zero crossings.

## `hopsize`

A parameter used to define the number audio of frames between frames size columns. 

# Example

```julia
using MusicProcessing

audio_one_channel = SampleBuf(rand(10000), 10)

zero_crossing_rate(audio_one_channel, 2048, 512)
```

"""
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
