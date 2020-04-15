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
include("polynomials/ChebyshevT.jl")

# compat; opt-in with `using Polynomials.PolyCompat`
include("polynomials/Poly.jl")
end # module
