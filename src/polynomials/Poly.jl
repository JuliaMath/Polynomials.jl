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
`polyint` to support older code. It should not be used for new code
base.

"""
struct Poly{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Poly(a::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where {T <: Number}
        # if a == [] we replace it with a = [0]
        Base.depwarn("Use of `Poly` from v1.0 forward will require `using Polynomials.PolyCompat`", :Poly)
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

Polynomials.domain(::Type{<:Poly}) = Polynomials.Interval(-Inf, Inf)
Polynomials.mapdomain(::Type{<:Poly}, x::AbstractArray) = x

function (p::Poly{T})(x::S) where {T,S}
    oS = one(x)
    length(p) == 0 && return zero(T) *  oS
    b = p[end]  *  oS
    @inbounds for i in (lastindex(p) - 1):-1:0
        b = p[i]*oS .+ x * b
    end
    return b
end


function Polynomials.fromroots(P::Type{<:Poly}, r::AbstractVector{T}; var::Polynomials.SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j in 1:n, i in j:-1:1
        c[i + 1] = c[i + 1] - r[j] * c[i]
    end
    return Poly(reverse(c), var)
end


function Polynomials.vander(P::Type{<:Poly}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    @inbounds for i in 1:n
        A[:, i + 1] = A[:, i] .* x
    end
    return A
end


function Polynomials.integrate(p::Poly{T}, k::S) where {T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(k)
        return Poly([NaN])
    end
    n = length(p)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return Poly(a2, p.var)
end


function Polynomials.derivative(p::Poly{T}, order::Integer = 1) where {T}
    order < 0 && error("Order of derivative must be non-negative")
    order == 0 && return p
    hasnan(p) && return Poly(T[NaN], p.var)
    order > length(p) && return zero(Poly{T})

    n = length(p)
    a2 = Vector{T}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return Poly(a2, p.var)
end


function Polynomials.companion(p::Poly{T}) where T
    d = length(p) - 1
    d < 1 && error("Series must have degree greater than 1")
    d == 1 && return diagm(0 => [-p[0] / p[1]])

    R = eltype(one(T) / p.coeffs[end])
    comp = diagm(-1 => ones(R, d - 1))
    monics = p.coeffs ./ p.coeffs[end]
    comp[:, end] .= -monics[1:d]
    return comp
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

function Base.divrem(num::Poly{T}, den::Poly{S}) where {T,S}
    num.var != den.var && error("Polynomials must have same variable")
    n = length(num) - 1
    m = length(den) - 1
    if m == 0 && den[0] ≈ 0 throw(DivideError()) end
    R = typeof(one(T) / one(S))
    P = Poly{R}
    deg = n - m + 1
    if deg ≤ 0
        return zero(P), convert(P, num)
    end
    q_coeff = zeros(R, deg)
    r_coeff = R.(num[0:n])
    @inbounds for i in n:-1:m
        q = r_coeff[i + 1] / den[m]
        q_coeff[i - m + 1] = q
        @inbounds for j in 0:m
            elem = den[j] * q
            r_coeff[i - m + j + 1] -= elem
        end
    end
    return P(q_coeff, num.var), P(r_coeff, num.var)
end

Polynomials.showterm(io::IO, ::Type{Poly{T}}, pj::T, var, j, first::Bool, mimetype) where {T} = showterm(io, Polynomial{T}, pj, var, j, first, mimetype)



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

# polyfit was deprecated to avoid a default calling `Poly`. Once
# PolyCompat is required, it can be used again
polyfit(x, y, n = length(x) - 1, sym=:x) = fit(Poly, x, y, n; var = sym)
polyfit(x, y, sym::Symbol) = fit(Poly, x, y, var = sym)

export Poly, poly, polyval, polyint, polyder#, polyfit


## Pade
include("../pade.jl")
using .PadeApproximation
export Pade
export padeval


end
