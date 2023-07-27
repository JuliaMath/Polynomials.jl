module Polynomials

#  using GenericLinearAlgebra ## remove for now. cf: https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl/pull/71#issuecomment-743928205
using LinearAlgebra
import Base: evalpoly
using Setfield

include("abstract.jl")
include("show.jl")
include("plots.jl")
include("contrib.jl")

# Interface for all AbstractPolynomials
include("common.jl")

# Polynomials
include("polynomials/standard-basis.jl")
include("polynomials/Polynomial.jl")
#include("polynomials/ImmutablePolynomial.jl")
#include("polynomials/SparsePolynomial.jl")
#include("polynomials/LaurentPolynomial.jl")

include("polynomials/pi_n_polynomial.jl")
include("polynomials/ngcd.jl")
include("polynomials/multroot.jl")
include("polynomials/factored_polynomial.jl")
include("polynomials/ChebyshevT.jl")

# polynomials with explicit basis
include("abstract-polynomial.jl")
include("basis-utils.jl")
include("polynomial-basetypes/mutable-dense-polynomial.jl")
include("polynomial-basetypes/immutable-dense-polynomial.jl")
include("polynomial-basetypes/mutable-sparse-polynomial.jl")



include("standard-basis/standard-basis.jl")
include("standard-basis/standard-dense.jl")
include("standard-basis/standard-immutable.jl")
include("standard-basis/standard-sparse.jl")

include("promotions.jl")

include("polynomials/chebyshev.jl") # wrong place


# Rational functions
include("rational-functions/common.jl")
include("rational-functions/rational-function.jl")
include("rational-functions/fit.jl")
#include("rational-functions/rational-transfer-function.jl")
include("rational-functions/plot-recipes.jl")

# compat; opt-in with `using Polynomials.PolyCompat`
include("polynomials/Poly.jl")

if !isdefined(Base, :get_extension)
    include("../ext/PolynomialsChainRulesCoreExt.jl")
    include("../ext/PolynomialsMakieCoreExt.jl")
    include("../ext/PolynomialsMutableArithmeticsExt.jl")
end

include("precompiles.jl")

end # module
