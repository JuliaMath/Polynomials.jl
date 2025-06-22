module PolynomialsMakieExt

using Polynomials
import Makie

function Makie.convert_arguments(P::Type{<:Makie.XYBased}, p::AbstractPolynomial)
    xs = Polynomials.poly_interval(p)
    return Makie.convert_arguments(P, xs, p.(xs))
end

Makie.plottype(p::AbstractPolynomial) = Makie.Lines

end
