import MakieCore

function MakieCore.convert_arguments(P::Type{<:MakieCore.XYBased}, p::AbstractPolynomial) 
    xs = poly_interval(p)
    return MakieCore.convert_arguments(P, xs, p.(xs))
end

MakieCore.plottype(p::AbstractPolynomial) = MakieCore.Lines
