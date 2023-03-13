module PolynomialsMakieCoreExt

using Polynomials
import MakieCore

function MakieCore.convert_arguments(P::Type{<:MakieCore.XYBased}, p::AbstractPolynomial)
    xs = Polynomials.poly_interval(p)
    return MakieCore.convert_arguments(P, xs, p.(xs))
end

MakieCore.plottype(p::AbstractPolynomial) = MakieCore.Lines

end
