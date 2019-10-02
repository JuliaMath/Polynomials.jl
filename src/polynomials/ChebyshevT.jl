

struct ChebyshevT{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
end