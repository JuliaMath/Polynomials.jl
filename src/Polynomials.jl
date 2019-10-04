module Polynomials

const SymbolLike = Union{AbstractString,Char,Symbol}
abstract type AbstractPolynomial{T <: Number} end

# Interface for all AbstractPolynomials
include("common.jl")

# New implementations
include("polynomials/Polynomial.jl")
include("polynomials/ChebyshevT.jl")
include("polynomials/ChebyshevU.jl")

# Deprecated -> Will be removed
include("polynomials/Poly.jl")

include("pade.jl")
include("deprecated.jl")

include("show.jl")
include("plots.jl")

end # module
