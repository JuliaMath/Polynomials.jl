

struct ChebyshevT{T <: Number} <: AbstractPolynomial
    coeffs::Vector{T}
    var::Symbol
end