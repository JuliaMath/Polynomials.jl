module PolyCompat

using ..Polynomials

#=
Compat support for old code. This will be opt-in by v1.0, through "using Polynomials.PolyCompat"
=#

## Poly

"""
    Poly{T}

Type of polynomial to support legacy code. Use of this type  is  not  encouraged.

This type provides support for `poly`, `polyval`, `polyder`, and
`polyint` to support older code. It should not be used for a new code
base. Call `using Polynomials.PolyCompat` to enable this module.

"""
struct Poly{T <: Number} <: Polynomials.StandardBasisPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Poly(a::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where {T <: Number}
        # if a == [] we replace it with a = [0]
        if length(a) == 0
            return new{T}(zeros(T, 1), Symbol(var))
        else
        # determine the last nonzero element and truncate a accordingly
            last_nz = findlast(!iszero, a)
            a_last = max(1, last_nz === nothing ? 0 : last_nz)
            new{T}(a[1:a_last], Symbol(var))
        end
    end
end

Polynomials.@register Poly

Base.convert(P::Type{<:Polynomial}, p::Poly{T}) where {T} = P(p.coeffs, p.var)

function (p::Poly{T})(x::S) where {T,S}
    oS = one(x)
    length(p) == 0 && return zero(T) *  oS
    b = p[end]  *  oS
    @inbounds for i in (lastindex(p) - 1):-1:0
        b = p[i]*oS .+ x * b
    end
    return b
end


function Base.:+(p1::Poly, p2::Poly)
    p1.var != p2.var && error("Polynomials must have same variable")
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n]
    return Poly(c, p1.var)
end

function Base.:*(p1::Poly{T}, p2::Poly{S}) where {T,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    n = length(p1) - 1
    m = length(p2) - 1
    R = promote_type(T, S)
    c = zeros(R, m + n + 1)
    for i in 0:n, j in 0:m
        c[i + j + 1] += p1[i] * p2[j]
    end
    return Poly(c, p1.var)
end


## Compat
## As this is an older package with many examples out in the wild
## rather than remove these, we limit them to this `Poly` type only

poly(r, var = :x) = fromroots(Poly, r; var = var)

polyval(p::Poly, x::Number) = p(x)
polyval(p::Poly, x) = p.(x)
polyval(p::AbstractPolynomial, x) = error("`polyval` is a legacy name for use with `Poly` objects only. Use `p(x)`.")

function Base.getproperty(p::Poly, nm::Symbol)
    if nm == :a
        return getfield(p, :coeffs)
    end
    return getfield(p, nm)
end

polyint(p::Poly, C = 0) = integrate(p, C)
polyint(p::Poly, a, b) = integrate(p, a, b)
polyint(p::AbstractPolynomial, args...)  = error("`polyint` is a legacy name for use with `Poly` objects only. Use `integrate(p,...)`.")

polyder(p::Poly, ord = 1) = derivative(p, ord)
polyder(p::AbstractPolynomial, args...) =  error("`polyder` is a legacy name for use with `Poly` objects only. Use `derivative(p,[order=1])`.")

polyfit(x, y, n = length(x) - 1, sym=:x) = fit(Poly, x, y, n; var = sym)
polyfit(x, y, sym::Symbol) = fit(Poly, x, y, var = sym)

export Poly, poly, polyval, polyint, polyder, polyfit

## Pade
include("../pade.jl")
using .PadeApproximation
export Pade
export padeval


end
