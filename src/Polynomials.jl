module Polynomials

using LinearAlgebra
using Intervals

include("abstract.jl")
include("show.jl")
include("plots.jl")
include("contrib.jl")

# Polynomials
include("polynomials/Polynomial.jl")

include("polynomials/ChebyshevT.jl")

# to be deprecated, then removed
include("polynomials/Poly.jl")
include("pade.jl")

# Interface for all AbstractPolynomials
include("common.jl")

end # module
