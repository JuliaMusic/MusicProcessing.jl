
export melspectrogram, mfcc

"""The Mel spectrogram"""
type MelSpectrogram{T, F <: Union{Frequencies,Range}} <: DSP.Periodograms.TFR{T}
    power::Matrix{T}
    mels::F
    time::FloatRange{Float64}
end
heatmap(tfr::MelSpectrogram) = log10(tfr.power)
DSP.freq(tfr::MelSpectrogram) = tfr.mels
Base.time(tfr::MelSpectrogram) = tfr.time

"""A TFR subtype where each column represents an MFCC vector"""
type MFCC{T, F <: Union{Frequencies,Range}} <: DSP.Periodograms.TFR{T}
    mfcc::Matrix{T}
    number::F
    time::FloatRange{Float64}
end
heatmap(tfr::MFCC) = tfr.mfcc
DSP.freq(tfr::MFCC) = tfr.number
Base.time(tfr::MFCC) = tfr.time

""""""
function hz_to_mel{F <: Real}(frequencies::Union{F, AbstractArray{F}})
    f_min = 0f0
    f_sp = 200f0 / 3

    mels = collect((frequencies - f_min) / f_sp)

    min_log_hz = 1000f0
    min_log_mel = (min_log_hz - f_min) / f_sp
    logstep = log(6.4f0) / 27f0

    @inbounds for i = 1:length(mels)
        if frequencies[i] >= min_log_hz
            mels[i] = min_log_mel + log(frequencies[i] / min_log_hz) / logstep
        end
    end

    mels
end

""""""
function mel_to_hz{F <: Real}(mels::Union{F, AbstractArray{F}})
    f_min = 0f0
    f_sp = 200f0 / 3
    frequencies = collect(f_min + f_sp * mels)

    min_log_hz = 1000f0
    min_log_mel = (min_log_hz - f_min) / f_sp
    logstep = log(6.4f0) / 27f0

    @inbounds for i = 1:length(frequencies)
        if mels[i] >= min_log_mel
            frequencies[i] = min_log_hz * exp(logstep * (mels[i] - min_log_mel))
        end
    end

    frequencies
end

""""""
function mel_frequencies(nmels::Int = 128, fmin::Real = 0.0f0, fmax::Real = 11025f0)
    min_mel = hz_to_mel(fmin)[1]
    max_mel = hz_to_mel(fmax)[1]

    mels = linspace(min_mel, max_mel, nmels)
    mel_to_hz(mels)
end

""""""
function mel(samplerate::Real, nfft::Int, nmels::Int = 128, fmin::Real = 0f0, fmax::Real = samplerate/2f0)
    weights = zeros(Float32, nmels, (nfft >> 1) + 1)
    fftfreqs = fft_frequencies(samplerate, nfft)
    melfreqs = mel_frequencies(nmels + 2, fmin, fmax)
    enorm = 2f0 ./ (melfreqs[3:end] - melfreqs[1:nmels])

    for i in 1:nmels
        lower = (fftfreqs - melfreqs[i]) / (melfreqs[i+1] - melfreqs[i])
        upper = (melfreqs[i+2] - fftfreqs) / (melfreqs[i+2] - melfreqs[i+1])

        weights[i, :] = max(0, min(lower, upper)) * enorm[i]
    end

    weights
end

function melspectrogram{T}(audio::SampleBuf{T, 1}, windowsize::Int = 1024, hopsize::Int = windowsize >> 2;
                        nmels::Int = 128, fmin::Real = 0f0, fmax::Real = audio.samplerate.val / 2f0)
    samplerate = audio.samplerate.val
    nfft = DSP.nextfastfft(windowsize)
    S = spectrogram(audio, windowsize, hopsize).power
    data = mel(samplerate, nfft, nmels, fmin, fmax) * S
    nframes = size(data, 2)
    MelSpectrogram(data, linspace(hz_to_mel(fmin)[1], hz_to_mel(fmax)[1], nmels), (0.0:nframes-1) * hopsize / samplerate)
end

function mfcc{T}(audio::SampleBuf{T, 1}, windowsize::Int = 1024, hopsize::Int = windowsize >> 2;
              nmfcc::Int = 20, nmels::Int = 128, fmin::Real = 0f0, fmax::Real = audio.samplerate.val / 2f0)
    if nmfcc >= nmels
        error("number of mfcc components should be less than the number of mel frequency bins")
    end

    M = melspectrogram(audio, windowsize, hopsize; nmels = nmels, fmin = fmin, fmax = fmax)
    mfcc = dct(nmfcc, nmels) * power(M)

    for frame in 1:size(mfcc, 2)
        mfcc[:, frame] /= norm(mfcc[:, frame])
    end

    MFCC(mfcc, 1:nmfcc, time(M))
end
