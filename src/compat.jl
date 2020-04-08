## We have renamed the MATLAB/numpy type names to more Julian names
## How to keep the old names during a transition is the question.
## The plan: keep these to ensure underlying changes are not disruptive
## For now  we ensure compatability by defining these for `Poly` objects such
## that they do not signal a deprecation (save polyfit)),
## but  do for other `AbstractPolynomial` types.
## At v1.0, it is likely these will be removed.


include("polynomials/Poly.jl")
using .PolyCompat
export Poly
export poly, polyval, polyint, polyder




## Pade
## Pade will  likely be moved into a separate pacakge
include("pade.jl")
using .PadeApproximation
export Pade
export padeval
