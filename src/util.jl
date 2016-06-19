# utility functions

hertz(value::Real) = value * Hz
hertz{F <: Real}(value::SIUnits.SIQuantity{F,0,0,-1,0,0,0,0,0,0}) = value

fft_frequencies(samplerate::Real, nfft::Int) = collect(linspace(0f0, samplerate / 2f0, (nfft >> 1) + 1))

""""""
function dct(nfilters::Int, ninput::Int)
    basis = Array(Float32, nfilters, ninput)
    samples = (1f0:2f0:2ninput) * π / 2ninput
    for i = 1:nfilters
        basis[i, :] = cos(i * samples)
    end

    basis *= sqrt(2f0/ninput)
    basis
end
