module PolynomialsFFTWExt

using Polynomials
import Polynomials: MutableDensePolynomial, StandardBasis, Pad
import FFTW
import FFTW: fft, ifft
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
