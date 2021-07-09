module PolyCompat

using ..Polynomials
indeterminate = Polynomials.indeterminate

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
struct Poly{T <: Number,X} <: Polynomials.StandardBasisPolynomial{T,X}
    coeffs::Vector{T}
    function Poly(a::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where {T <: Number}
        # if a == [] we replace it with a = [0]
        X = Symbol(var)
        if length(a) == 0
            return new{T,X}(zeros(T, 1))
        else
        # determine the last nonzero element and truncate a accordingly
            last_nz = findlast(!iszero, a)
            a_last = max(1, last_nz === nothing ? 0 : last_nz)
            new{T,X}(a[1:a_last])
        end
    end
end

Polynomials.@register Poly

Poly{T,X}(coeffs::AbstractVector{S}) where {T,X,S} = Poly(convert(Vector{T},coeffs), X)

Base.convert(P::Type{<:Polynomial}, p::Poly{T}) where {T} = Polynomial{eltype(P),indeterminate(P,p)}(p.coeffs)
Base.convert(P::Type{<:Poly}, p::Poly{T,X}) where {T<:Number,X} = Polynomials.constructorof(P){_eltype(P),indeterminate(P,p)}(p.coeffs)
Base.eltype(P::Type{<:Poly{T,X}}) where {T, X} = P
_eltype(::Type{<:Poly{T}}) where  {T} = T
_eltype(::Type{Poly}) =  Float64

# when interating over poly return monomials
function Base.iterate(p::Poly, state=nothing)
    i = 0
    state == nothing && return (p[i]*one(p), i)
    j = degree(p)
    s = state + 1
    i <= state < j && return (p[s]*Polynomials.basis(p,s), s)
    return nothing
end
Base.collect(p::Poly) = [pᵢ for pᵢ ∈ p]

# need two here as `eltype(P)` is `_eltype(P)`.
Base.zero(::Type{P}) where {P <: Poly} = Poly{_eltype(P), Polynomials.indeterminate(P)}([0])
Base.zero(::Type{P},var::Polynomials.SymbolLike) where {P <: Poly} = Poly(zeros(_eltype(P),1), var)
Base.one(::Type{P}) where {P <: Poly} = Poly{_eltype(P), Polynomials.indeterminate(P)}([1])
Base.one(::Type{P},var::Polynomials.SymbolLike) where {P <: Poly} = Poly(ones(_eltype(P),1), var)
Polynomials.variable(::Type{P}) where {P <: Poly} = Poly{_eltype(P), Polynomials.indeterminate(P)}([0,1])
Polynomials.variable(::Type{P},var::Polynomials.SymbolLike) where {P <: Poly} = Poly(_eltype(P)[0,1], var)
function Polynomials.basis(P::Type{<:Poly}, k::Int, _var::Polynomials.SymbolLike=:x; var=_var)
    zs = zeros(Int, k+1)
    zs[end] = 1
    Polynomials.constructorof(P){_eltype(P), Symbol(var)}(zs)
end

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
    indeterminate(p1) != indeterminate(p2) && error("Polynomials must have same variable")
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n-1]
    return Poly(c, indeterminate(p1))
end

function Base.:*(p1::Poly{T}, p2::Poly{S}) where {T,S}
    indeterminate(p1) != indeterminate(p2) && error("Polynomials must have same variable")
    n = length(p1) - 1
    m = length(p2) - 1
    R = promote_type(T, S)
    c = zeros(R, m + n + 1)
    for i in 0:n, j in 0:m
        c[i + j + 1] += p1[i] * p2[j]
    end
    return Poly(c, indeterminate(p1))
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

function Polynomials.integrate(p::P, k::S) where {T, X, P <: Poly{T, X}, S<:Number}

    R = eltype((one(T)+one(S))/1)
    Q = Poly{R,X}
    if hasnan(p) || isnan(k)
        return P([NaN]) # keep for Poly, not Q
    end
    n = length(p)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return Q(a2)
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
