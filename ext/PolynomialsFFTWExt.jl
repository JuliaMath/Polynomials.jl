module PolynomialsFFTWExt

using Polynomials
import Polynomials: MutableDensePolynomial, StandardBasis
import FFTW
import FFTW: fft, ifft

struct Pad{T} <: AbstractVector{T}
    a::Vector{T}
    n::Int
end
Base.length(a::Pad) = a.n
Base.size(a::Pad) = (a.n,)
function Base.getindex(a::Pad, i)
    u = length(a.a)
    i ≤ u && return a.a[i]
    return zero(first(a.a))
end

function Polynomials.poly_multiplication_fft(p::P, q::Q) where {B <: StandardBasis,X,
                                                                T <: AbstractFloat, P<:MutableDensePolynomial{B,T,X},
                                                                S <: AbstractFloat, Q<:MutableDensePolynomial{B,S,X}}

    N = 1 + degree(p) + degree(q)
    as = Pad(p.coeffs, N)
    bs = Pad(q.coeffs, N)
    us = fft(as)
    vs = fft(bs)
    cs = ifft(us .* vs)
    map!(real, cs, cs)
    MutableDensePolynomial{B, eltype(cs), X}(cs)

end



end
