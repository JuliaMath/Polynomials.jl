struct ChebyshevU{T <: Number} <: AbstractPolynomial
    coeffs::Vector{T}
    var::Symbol
end