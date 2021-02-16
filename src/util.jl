# utility functions

tofloat(array::AbstractArray{T}) where T = map(Float32, array)
tofloat(array::AbstractArray{Float32}) = array

fft_frequencies(samplerate::Real, nfft::Int) = collect(range(0f0, samplerate / 2f0, length=(nfft >> 1) + 1))


"""returns the number of frames when the signal is partitioned into overlapping frames"""
nframes(length::Int, framesize::Int, hopsize::Int) = div(length - framesize, hopsize) + 1

"""Provides a view of a signal partitioned into overlapping frames, similar to MATLAB's buffer"""
struct FrameView{T<:AbstractVector} <: AbstractVector{Vector}
    original::T
    framesize::Int
    hopsize::Int
    nframes::Int

    function FrameView(original, framesize, hopsize)
        # n = noverlap is a problem - the algorithm will not terminate.
        #@boundscheck((0 < hopsize <= framesize), error("hopsize must be between zero and framesize"))
        new{T}(original, framesize, hopsize, nframes(length(original), framesize, hopsize))
    end
end

function Base.getindex(x::FrameView, i::Int)
    #@boundscheck((i >= 1 && i <= x.k), Base.throw_boundserror(x, i))
    slice(x.original, (i-1) * x.hopsize + (1:x.framesize))
end

Base.size(x::FrameView) = (x.nframes,)

function Base.iterate(x::FrameView, i::Int)
   if i > x.nframes
       nothing
   else
       (x[i], i+1)
   end
end



"""returns the DCT filters"""
function dct(nfilters::Int, ninput::Int)
    basis = Array{Float32}(undef, nfilters, ninput)
    samples = (1f0:2f0:2ninput) * π / 2ninput
    for i = 1:nfilters
        for j in 1:length(samples)
            basis[i, j] = cos(i * samples[j])
        end
    end

    basis *= sqrt(2f0/ninput)
    basis
end

"""returns the number of times that the value changes its sign"""
function zero_crossings(array::AbstractVector{T}, size::Int = length(array), offset::Int = 0) where {T<:Real}
    result = 0
    previous = 0
    for i in offset .+ (1:size)
        number = array[i]
        sgn = number == 0 ? 0 : number > 0 ? 1 : -1
        if sgn != previous
            result += 1
        end
        previous = sgn
    end
    result
end
