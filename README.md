# MusicProcessing.jl

<!--[![Build Status](https://travis-ci.org/jongwook/MusicProcessing.jl.svg?branch=master)](https://travis-ci.org/jongwook/MusicProcessing.jl)-->

MusicProcessing.jl is a music and audio processing library for Julia, inspired by [librosa](github.com/librosa/librosa). It is not feature complete and in a very early stage of development.

# Performance

Thanks to Julia's performance optimizations, it is significantly faster than librosa, a mature library in Python

![Imgur](http://i.imgur.com/YHoALHr.png)

All measurements are done by averaging over 100 repetitions, after one warmup run.

# Usage

The following commands will display a graphic visualization and/or an HTML5 `<audio>` component for playing audio, when run in IJulia.

### Loading an audio file

```julia
julia> using MusicProcessing, FileIO, MP3
julia> audio = load("Sour_Tennessee_Red_Sting.mp3")
```

### Converting to a mono audio

```julia
julia> audio = mono(audio)
```

### Resampling in 22050 Hz

```julia
julia> audio = resample(audio, 22050Hz)
```

### Speeding up the audio

```julia
julia> speedup(audio, 2)
```

### Pitch-shifting

```julia
julia> pitchshift(audio, 4)
```

### Displaying Spectrogram

```julia
julia> spectrogram(audio)
```

### Displaying Mel Spectrogram

```julia
julia> melspectrogram(audio)
```

### Displaying MFCC

```julia
julia> mfcc(audio)
```

# Roadmap

There are a lot to be implemented, including and not limited to:

* Harmonic Features
 * CQT, Chroma, Tonnetz
* Rhythmic Features
 * Onset Detection, Beat Detection
* Melodic Features
 * F0 tracking, multi-pitch tracking
* Source Separation
 * Harmonic-Percussive Source Separation
* Performance Tuning, Tests...
