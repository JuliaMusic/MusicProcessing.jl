
import Base.writemime
export writemime

# methods to translate

"""return the frequency ticks of a spectrogram, rounded to the nearest integers"""
function yticklabels(tfr::DSP.Periodograms.Spectrogram, yticks::Array)
    map(Int, map(round, yticks)), "frequency (Hz)"
end

"""return Hz values corresponding to given mel frequencies"""
function yticklabels(tfr::MelSpectrogram, yticks::Array)
    map(Int, map(round, mel_to_hz(yticks))), "frequency (Hz)"
end

"""return MFCC coefficient numbers"""
function yticklabels(tfr::MFCC, yticks::Array)
    map(Int, map(round, yticks)), "MFCC number"
end

# writemime is called by

"""Display the waveform of given audio, along with an HTML audio element with controls"""
function writemime{T, N}(io::IO, mime::MIME"text/html", x::SampleBuf{T, N, Hertz})
    nchannels = SampledSignals.nchannels(x)
    nframes = SampledSignals.nframes(x)
    samplerate = SampledSignals.samplerate(x).val

    write(io, """<p>
        $(@sprintf("%.2f", duration(x).val))s
        $((nchannels == 2) ? "stereo" : (nchannels == 1) ? "mono" : "$nchannels-channel") audio;
        sample rate = $samplerate Hz,
    </p>""")

    @eval import MP3, PyPlot

    PyPlot.ioff()
    PyPlot.figure(figsize = (8, 2 * nchannels))

    for i = 1:nchannels
        PyPlot.subplot(nchannels, 1, i)
        PyPlot.plot((1.0:nframes)/samplerate, x.data[:, i])
        PyPlot.gca()[:spines]["top"][:set_visible](false)
        PyPlot.gca()[:spines]["left"][:set_visible](false)
        PyPlot.gca()[:spines]["right"][:set_visible](false)
        if i != nchannels
            PyPlot.gca()[:get_xaxis]()[:set_visible](false)
            PyPlot.gca()[:spines]["bottom"][:set_visible](false)
        else
            PyPlot.xlabel("time (seconds)")
        end
        PyPlot.gca()[:get_yaxis]()[:set_visible](false)
    end

    write(io, """
        <p>
            <img src="data:image/png;base64,""")

    base64 = Base64EncodePipe(io)
    writemime(base64, MIME("image/png"), PyPlot.gcf())
    close(base64)

    write(io, """">
        </p>
    """)

    lame = MP3.lame_init()
    try
        info = MP3.MP3INFO(x)
        MP3.lame_set_num_samples(lame, info.nframes)
        MP3.lame_set_in_samplerate(lame, info.samplerate)
        if info.nchannels == 1
            MP3.lame_set_num_channels(lame, 1)
            MP3.lame_set_mode(lame, MP3.LAME_MONO)
        else
            MP3.lame_set_num_channels(lame, 2)
            MP3.lame_set_mode(lame, MP3.LAME_JOINT_STEREO)
        end
        MP3.lame_set_out_samplerate(lame, info.samplerate)
        MP3.lame_set_brate(lame, 320)
        MP3.lame_init_params(lame)

        base64 = Base64EncodePipe(io)
        sink = MP3.MP3FileSink(lame, info, base64)

        write(io, """
            <audio controls="controls" style="width:800px;">
                <source type="audio/mp3" src="data:audio/mp3;base64,""")

        @eval import MP3.unsafe_write
        write(sink, x)
        close(base64)

        write(io, """"
                Your browser does not support the audio element.
            </audio>
        """)
    finally
        MP3.lame_close(lame)
    end
end

heatmap(tfr::DSP.Periodograms.TFR) = log10(power(tfr))

function draw_heatmap(tfr::DSP.Periodograms.TFR)
    X = time(tfr)
    Y = freq(tfr)
    Z = heatmap(tfr)

    PyPlot.pcolormesh(X, Y, Z)
    PyPlot.xlim(first(X), last(X))
    PyPlot.ylim(first(Y), last(Y))

    (yticks, ylabel) = yticklabels(tfr, PyPlot.yticks()[1])
    PyPlot.gca()[:set_yticklabels](yticks)
    PyPlot.gca()[:spines]["top"][:set_visible](false)
    PyPlot.gca()[:spines]["right"][:set_visible](false)

    PyPlot.xlabel("time (seconds)")
    PyPlot.ylabel(ylabel)
end

"""Display a spectrogram"""
function writemime{R <: DSP.Periodograms.TFR}(io::IO, mime::MIME"image/png", tfr::R)
    @eval import PyPlot
    PyPlot.ioff()
    PyPlot.figure(figsize=(8, 4))
    draw_heatmap(tfr)
    writemime(io, mime, PyPlot.gcf())
end

"""Display multichannel spectrogram"""
function writemime{R <: DSP.Periodograms.TFR}(io::IO, mime::MIME"image/png", tfrs::Array{R, 1})
    nchannels = length(tfrs)

    @eval import PyPlot
    PyPlot.ioff()
    PyPlot.figure(figsize=(8, 8))

    for i = 1:nchannels
        PyPlot.subplot(nchannels, 1, i)
        draw_heatmap(tfrs[i])

        if i != nchannels
            PyPlot.gca()[:get_xaxis]()[:set_visible](false)
            PyPlot.gca()[:spines]["bottom"][:set_visible](false)
        end
    end
    writemime(io, mime, PyPlot.gcf())
end
