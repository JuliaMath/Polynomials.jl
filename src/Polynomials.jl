module Polynomials

using LinearAlgebra
using Intervals

include("abstract.jl")
include("show.jl")
include("plots.jl")

# Polynomials
include("polynomials/Polynomial.jl")
include("polynomials/ChebyshevT.jl")
include("polynomials/ChebyshevU.jl")

include("polynomials/Poly.jl") # Deprecated -> Will be removed
include("pade.jl")
include("deprecated.jl")


# Interface for all AbstractPolynomials
include("common.jl")

end # module
