module Polynomials

#  using GenericLinearAlgebra ## remove for now. cf: https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl/pull/71#issuecomment-743928205
using LinearAlgebra
using Intervals

if VERSION >= v"1.4.0"
    import Base: evalpoly
end

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
include("polynomials/pi_n_polynomial.jl")
include("polynomials/ngcd.jl")
include("polynomials/multroot.jl")

include("polynomials/ChebyshevT.jl")

# Rational functions
if VERSION >= v"1.2.0"
    include("rational-functions/common.jl")
    include("rational-functions/rational-function.jl")
    include("rational-functions/fit.jl")
    #include("rational-transfer-function.jl")
    include("rational-functions/plot-recipes.jl")
end


# compat; opt-in with `using Polynomials.PolyCompat`
include("polynomials/Poly.jl")

end # module
