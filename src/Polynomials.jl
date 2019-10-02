module Polynomials

# Interface for all AbstractPolynomials
include("common.jl")

# Original (to deprecate) code
include("old.jl")
include("pade.jl")
include("show.jl")
include("PlotRecipes.jl")

# New implementations
include("polynomials/Polynomial.jl")
include("polynomials/ChebyshevT.jl")
include("polynomials/ChebyshevU.jl")

end # module
