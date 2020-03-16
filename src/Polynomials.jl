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
#include("polynomials/ChebyshevU.jl")


include("polynomials/Poly.jl") # to be deprecated, then removed
include("pade.jl")
include("compat.jl") # Where we keep deprecations

# Interface for all AbstractPolynomials
include("common.jl")

end # module
