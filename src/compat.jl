## We have renamed the MATLAB/numpy type names to more Julian names
## How to keep the old names during a transition is the question.
## The plan: keep these to ensure underlying changes are not disruptive
## For now  we ensure compatability by defining these for `Poly` objects such
## that they do not signal a deprecation (save polyfit)),
## but  do for other `AbstractPolynomial` types.
## At v1.0, it is likely these will be removed.

## Ensure compatability for now
@deprecate polyval(p::AbstractPolynomial, x::Number)  p(x)
@deprecate polyval(p::AbstractPolynomial, x)  p.(x)


@deprecate polyint(p::AbstractPolynomial, C = 0)  integrate(p, C)
@deprecate polyint(p::AbstractPolynomial, a, b)  integrate(p, a, b)

@deprecate polyder(p::AbstractPolynomial, ord = 1)  derivative(p, ord)

@deprecate polyfit(x, y, n = length(x) - 1, sym=:x)  fit(Poly, x, y, n; var = sym)
@deprecate polyfit(x, y, sym::Symbol)  fit(Poly, x, y, var = sym)


include("polynomials/Poly.jl")
using .PolyCompat
export Poly
export poly, polyval, polyint, polyder, polyfit




## Pade
## Pade will  likely be moved into a separate pacakge
include("pade.jl")
using .PadeApproximation
export Pade
export padeval
