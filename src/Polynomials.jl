module Polynomials

using LinearAlgebra
using Intervals
using OffsetArrays

include("abstract.jl")
include("show.jl")
include("plots.jl")
include("contrib.jl")

# Interface for all AbstractPolynomials
include("common.jl")


# Polynomials
include("polynomials/standard-basis.jl")
include("polynomials/Polynomial.jl")
include("polynomials/ImmutablePolynomial.jl")
include("polynomials/SparsePolynomial.jl")
include("polynomials/LaurentPolynomial.jl")

include("polynomials/ChebyshevT.jl")

# compat; opt-in with `using Polynomials.PolyCompat`
include("polynomials/Poly.jl")

end # module
