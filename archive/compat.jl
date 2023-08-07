## We have renamed the MATLAB/numpy type names to more Julian names
## How to keep the old names during a transition is the question.
## The plan: keep these to ensure underlying changes are not disruptive
## For now  we ensure compatibility by defining these for `Poly` objects such
## that they do not signal a deprecation (save polyfit)),
## but  do for other `AbstractPolynomial` types.
## At v1.0, these will be opt-in via `using Polynomials.PolyCompat`


include("polynomials/Poly.jl")
using .PolyCompat
export Poly
export poly, polyval, polyint, polyder

## Pade
## Pade will  likely be moved into a separate package, for now we will put into PolyCompat
export Pade, padeval
