module Polynomials

const SymbolLike = Union{AbstractString,Char,Symbol}
abstract type AbstractPolynomial{T <: Number} end


include("common.jl")

# New implementations
include("polynomials/Polynomial.jl")
include("polynomials/ChebyshevT.jl")
include("polynomials/ChebyshevU.jl")

# Interface for all AbstractPolynomials
include("deprecated.jl")

# Original (to deprecate) code
# include("old.jl")
# include("pade.jl")
# include("show.jl")
# include("PlotRecipes.jl")

end # module
