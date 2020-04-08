module Polynomials

using LinearAlgebra
using Intervals

include("abstract.jl")
include("show.jl")
include("plots.jl")
include("contrib.jl")

# Interface for all AbstractPolynomials
include("common.jl")


# Polynomials
include("polynomials/Polynomial.jl")
include("polynomials/ImmutablePolynomial.jl")
include("polynomials/ChebyshevT.jl")

# to be deprecated, then removed
include("compat.jl")

end # module
