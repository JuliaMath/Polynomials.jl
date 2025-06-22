"""
    Polynomials.jl

Basic arithmetic, integration, differentiation, evaluation, root finding, and fitting for [univariate polynomials](https://en.wikipedia.org/wiki/Polynomial) in [Julia](https://julialang.org/).
"""
module Polynomials

#  using GenericLinearAlgebra ## remove for now. cf: https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl/pull/71#issuecomment-743928205
using LinearAlgebra
import Base: evalpoly
using Setfield
using SparseArrays
using OrderedCollections

include("abstract.jl")
include("show.jl")
include("plots.jl")
include("contrib.jl")

# Interface for all AbstractPolynomials
include("common.jl")

# polynomials with explicit basis
include("abstract-polynomial.jl")
include("polynomial-container-types/mutable-dense-polynomial.jl")
include("polynomial-container-types/mutable-dense-view-polynomial.jl")
include("polynomial-container-types/mutable-dense-laurent-polynomial.jl")
include("polynomial-container-types/immutable-dense-polynomial.jl")
include("polynomial-container-types/mutable-sparse-polynomial.jl")
include("polynomial-container-types/mutable-sparse-vector-polynomial.jl")
const PolynomialContainerTypes = (:MutableDensePolynomial, :MutableDenseViewPolynomial, :ImmutableDensePolynomial,
                                  :MutableDenseLaurentPolynomial, :MutableSparsePolynomial, :MutableSparseVectorPolynomial) # useful for some purposes
const ZeroBasedDensePolynomialContainerTypes = (:MutableDensePolynomial, :MutableDenseViewPolynomial, :ImmutableDensePolynomial)

include("polynomials/standard-basis/standard-basis.jl")
include("polynomials/standard-basis/polynomial.jl")
include("polynomials/standard-basis/pn-polynomial.jl")
include("polynomials/standard-basis/laurent-polynomial.jl")
include("polynomials/standard-basis/immutable-polynomial.jl")
include("polynomials/standard-basis/sparse-polynomial.jl")
include("polynomials/standard-basis/sparse-vector-polynomial.jl")

include("polynomials/ngcd.jl")
include("polynomials/multroot.jl")

include("polynomials/factored_polynomial.jl")
include("polynomials/chebyshev.jl")

include("promotions.jl")



# Rational functions
include("rational-functions/common.jl")
include("rational-functions/rational-function.jl")
include("rational-functions/fit.jl")
#include("rational-functions/rational-transfer-function.jl")
include("rational-functions/plot-recipes.jl")

# compat; opt-in with `using Polynomials.PolyCompat`
include("legacy/misc.jl")
include("legacy/Poly.jl")

include("precompiles.jl")

@static if !isdefined(Base, :get_extension)
    using Requires
end

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/PolynomialsCoreExt.jl")
    end
end

end # module
