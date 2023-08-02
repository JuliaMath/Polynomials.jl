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
include("polynomials/factored_polynomial.jl")

# polynomials with explicit basis
include("abstract-polynomial.jl")
include("polynomial-basetypes/mutable-dense-polynomial.jl")
include("polynomial-basetypes/mutable-dense-view-polynomial.jl")
include("polynomial-basetypes/mutable-dense-laurent-polynomial.jl")
include("polynomial-basetypes/immutable-dense-polynomial.jl")
include("polynomial-basetypes/mutable-sparse-polynomial.jl")

include("standard-basis/standard-basis.jl")
include("standard-basis/polynomial.jl")
include("standard-basis/pn-polynomial.jl")
include("standard-basis/laurent-polynomial.jl")
include("standard-basis/immutable-polynomial.jl")
include("standard-basis/sparse-polynomial.jl")

include("promotions.jl")

include("polynomials/ngcd.jl")
include("polynomials/multroot.jl")

include("polynomials/chebyshev.jl")


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
