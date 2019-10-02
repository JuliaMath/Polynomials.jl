# Poly type manipulations

module Polynomials
#todo: sparse polynomials?

include("common.jl")
include("polynomials/old.jl")

include("polynomials/Polynomial.jl")
include("polynomials/ChebyshevT.jl")

### Pull in others
include("show.jl") # display polynomials.
include("pade.jl")
include("PlotRecipes.jl")

end # module Poly
