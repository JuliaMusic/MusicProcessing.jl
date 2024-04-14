# complex number utils
import FFTW.libfftw3, FFTW.libfftw3f

function cis!(dest::AbstractVector{Complex{T}}, v::AbstractVector{T}) where {T<:Number}
    @inbounds for i in 1:length(v)
        dest[i] = cis(v[i])
    end
end

function multiply!(dest::AbstractVector{Complex{T}}, left::AbstractVector{Complex{T}}, right::AbstractVector{T}) where {T<:Number}
    @inbounds for i in 1:length(dest)
        dest[i] = left[i] * right[i]
    end
end

const ArrayLike{T} = Union{Vector{T}, Ptr{T}}

function plan_dft_c2r_1d!(dest::ArrayLike{Float32}, src::ArrayLike{Complex{Float32}}, nfft::Int)
    ccall((:fftwf_plan_dft_c2r_1d, libfftw3f.x), FFTW.PlanPtr, (Cint, Ptr{Complex{Float32}}, Ptr{Float32}, Cuint), nfft, src, dest, 0)
end

function plan_dft_c2r_1d!(dest::ArrayLike{Float64}, src::ArrayLike{Complex{Float64}}, nfft::Int)
    ccall((:fftw_plan_dft_c2r_1d, libfftw3.x), FFTW.PlanPtr, (Cint, Ptr{Complex{Float64}}, Ptr{Float64}, Cuint), nfft, src, dest, 0)
end

function execute_dft_c2r!(dest::ArrayLike{Float32}, src::ArrayLike{Complex{Float32}}, plan::FFTW.PlanPtr)
    ccall((:fftwf_execute_dft_c2r, libfftw3f.x), Nothing, (FFTW.PlanPtr, Ptr{Complex{Float32}}, Ptr{Float32}), plan, src, dest)
end

function execute_dft_c2r!(dest::ArrayLike{Float64}, src::ArrayLike{Complex{Float64}}, plan::FFTW.PlanPtr)
    ccall((:fftw_execute_dft_c2r, libfftw3.x), Nothing, (FFTW.PlanPtr, Ptr{Complex{Float64}}, Ptr{Float64}), plan, src, dest)
end

"""Performs inverse real-data FFT without allocating the result array"""
function irfft!(dest::Array{T}, src::ArrayLike{Complex{T}}, nfft::Int) where {T<:Number}
    plan = plan_dft_c2r_1d!(dest, src, nfft)
    execute_dft_c2r!(dest, src, plan)
    for i = 1:nfft
        dest[i] /= nfft
    end
end
