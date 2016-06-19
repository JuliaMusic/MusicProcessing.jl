
import Base.eps
import DSP.spectrogram, DSP.stft
export spectrogram, stft, istft, phase_vocoder

""""""
function spectrogram{T}(audio::SampleBuf{T, 1, Hertz},
                        windowsize::Int = 1024,
                        hopsize::Int = windowsize >> 2;
                        window = hanning, kwargs...)
    noverlap = windowsize - hopsize
    data = map(Float32, audio.data)
    DSP.spectrogram(data, windowsize, noverlap; fs = audio.samplerate.val, window = window, kwargs...)
end

""""""
function spectrogram{T}(audio::SampleBuf{T, 2, Hertz},
                        windowsize::Int = 1024,
                        hopsize::Int = windowsize >> 2;
                        window = hanning, kwargs...)
    noverlap = windowsize - hopsize
    data = map(Float32, audio.data)
    (mapslices(data, 1) do data
        DSP.spectrogram(data, windowsize, noverlap; fs = audio.samplerate.val, window = window, kwargs...)
    end)[:]
end

""""""
function spectrogram{T, N}(audio::SampleBuf{T, N, Hertz},
                           windowsize::Seconds,
                           hopsize::Seconds = windowsize / 4;
                           kwargs...)
    w = round(Int, windowsize * audio.samplerate)
    h = round(Int, hopsize * audio.samplerate)
    spectrogram(audio, w, h; kwargs...)
end


""""""
function stft{T}(audio::SampleBuf{T, 1, Hertz},
                 windowsize::Int = 1024,
                 hopsize::Int = windowsize >> 2;
                 window = hanning, kwargs...)
    noverlap = windowsize - hopsize
    data = map(Float32, audio.data)
    DSP.stft(data, windowsize, noverlap; window = window, kwargs...)
end

""""""
function stft{T}(audio::SampleBuf{T, 2, Hertz},
                 windowsize::Int = 1024,
                 hopsize::Int = windowsize >> 2; kwargs...)
    nchannels = SampledSignals.nchannels(audio)
    noverlap = windowsize - hopsize

    stft = Array(Matrix, nchannels)
        data = map(Float32, audio.data)
    for i in 1:nchannels
        stft[i] = DSP.stft(data[:, i], windowsize, noverlap; kwargs...)
    end
    cat(3, stft...)
end

""""""
function stft{T, N}(audio::SampleBuf{T, N, Hertz},
                    windowsize::Seconds,
                    hopsize::Seconds = windowsize / 4; kwargs...)
    w = round(Int, windowsize * audio.samplerate)
    h = round(Int, hopsize * audio.samplerate)
    stft(audio, w, h; kwargs...)
end


# functions to make window vector from a window function, vector or nothing.
@inline w(window::Function, size) = window(size)
@inline function w(window::AbstractVector, size)
    if size != length(window)
        error("the length of window does not match windowsize")
    end
    window
end
@inline w(window::Void,size) = 1


"""
    istft(stft, samplerate, windowsize, hopsize; window)

Reconstruct a mono audio from its STFT.

# Arguments
* `stft::Array{Complex{T}, 2}`: the STFT matrix, usually having (1+nfft/2) rows
* `samplerate::Union{Real, Hertz}`: Sample rate of the resulting audio
* `windowsize::Int`: window size used for STFT. (default: nfft)
* `hopsize::Int`: number of frames between STFT columns (default: windowsize/4)
* `nfft::Int`: the FFT size (default: 2*(size(stft,1)-1))
* `window::Union{Function, AbstractVector, Void}`: function or vector that has the window used for STFT (default: uniform window)
"""
function istft{T <: AbstractFloat}(stft::Array{Complex{T}, 2},
                                   samplerate::Union{Real, Hertz},
                                   windowsize::Int = 2 * (size(stft, 1) - 1),
                                   hopsize::Int = windowsize >> 2;
                                   nfft::Int = windowsize,
                                   window::Union{Function, AbstractVector, Void} = hanning)

    if windowsize > nfft
      error("window size should be less than or equal to nfft")
    end

    window = w(window, windowsize)

    columns = size(stft, 2)
    audio = zeros(T, windowsize + hopsize * (columns - 1))
    weights = zeros(audio)

    for i in 1:columns
        range = (i-1) * hopsize + (1:windowsize)
        spectrum = irfft(stft[:, i], nfft)[1:windowsize]

        audio[range] += window .* spectrum
        weights[range] += window
    end

    @inbounds for i in 1:length(audio)
        if weights[i] > eps(T)
            audio[i] /= weights[i]
        end
    end

    SampleBuf{T, 1, Hertz}(audio, hertz(samplerate))
end

"""
    istft(stft, samplerate, windowsize, hopsize; window)

Reconstruct a multichannel audio from its STFT.

# Arguments
* `stft::Array{Complex{T}, 3}`: the STFT matrix, usually having (1+nfft/2) rows
* `samplerate::Union{Real, Hertz}`: Sample rate of the resulting audio
* `windowsize::Int`: window size used for STFT. (default: nfft)
* `hopsize::Int`: number of frames between STFT columns (default: windowsize/4)
* `nfft::Int`: the FFT size (default: 2*(size(stft,1)-1))
* `window::Union{Function, AbstractVector, Void}`: function or vector that has the window used for STFT (default: uniform window)
"""
function istft{T <: AbstractFloat}(stft::Array{Complex{T}, 3},
                                   samplerate::Union{Real, Hertz},
                                   windowsize::Int = 2 * (size(stft, 1) - 1),
                                   hopsize::Int = windowsize >> 2; kwargs...)
    nchannels = size(stft, 3)
    buffers = cell(nchannels)

    # run ISTFT for each channel
    for i = 1:nchannels
        buffers[i] = istft(stft[:, :, i], samplerate, windowsize, hopsize; kwargs...)
    end

    nsamples = size(buffers[1], 1)
    audio = Array(T, nsamples, nchannels)

    for i = 1:nchannels
        audio[:, i] = buffers[i].data
    end

    # concatenate channels
    SampleBuf{T, 2, Hertz}(
        audio,
        hertz(samplerate)
    )
end


"""
    phase_vocoder(stft, speed, hopsize)

Phase vocoder. Given an STFT matrix, speed it up by a factor.
"""
function phase_vocoder{T <: AbstractFloat}(stft::Array{Complex{T}, 2},
                                           speed::Real,
                                           hopsize::Int = (size(stft, 1) - 1) >> 1)

    nbins = size(stft, 1)
    nframes = size(stft, 2)
    nfft = (nbins - 1) * 2
    timesteps = 1.0:speed:nframes
    stretched = zeros(Complex{T}, nbins, length(timesteps))
    phase_advance = linspace(0, π * hopsize, nbins)
    phase_acc = angle(stft[:, 1])

    @inbounds for (t, step) in enumerate(timesteps)
        left = floor(Int, step)
        right = left + 1
        columns = stft[:, left:min(right, nframes)]

        # weighting for linear magnitude interpolation
        alpha = step % 1.0
        if size(columns, 2) == 2
            mag = ((1.0 - alpha) * abs(columns[:, 1]) + alpha * abs(columns[:, 2]))
        else
            mag = columns[:, 1]
        end

        # store to output array
        stretched[:, t] = mag .* exp(phase_acc * im)

        # nothing to do more
        if t == size(stretched, 2)
            break
        end

        # compute phase advance
        dphase = angle(columns[:, 2]) - angle(columns[:, 1]) - phase_advance

        # wrap to -pi:pi range
        dphase -= 2π * round(dphase / 2π)

        # accumulate phase
        phase_acc += phase_advance + dphase
    end

    stretched
end

"""
    phase_vocoder(stft, speed, hopsize)

Phase vocoder. Given an STFT matrix, speed it up by a factor.
"""
function phase_vocoder{T <: AbstractFloat}(stft::Array{Complex{T}, 3},
                                           speed::Real,
                                           hopsize::Int = (size(stft, 1) - 1) >> 1)
    mapslices(stft, 1:2) do stft
        phase_vocoder(stft, speed, hopsize)
    end
end
