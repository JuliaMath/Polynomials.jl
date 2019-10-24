export Polynomial

"""
    Polynomial{T<:Number}(coeffs::AbstractVector{T}, var=:x)

Construct a polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct this through
`Polynomial([a_0, a_1, ..., a_n])`.

The usual arithmetic operators are overloaded to work with polynomials as well as
with combinations of polynomials and scalars. However, operations involving two
polynomials of different variables causes an error.

# Examples

```jldoctest
julia> Polynomial([1, 0, 3, 4])
Polynomial(1 + 3*x^2 + 4*x^3)

julia> Polynomial([1, 2, 3], :s)
Polynomial(1 + 2*s + 3*s^2)

julia> one(Polynomial)
Polynomial(1.0)
```
"""
struct Polynomial{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Polynomial{T}(coeffs::Vector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

@register Polynomial

domain(::Type{<:Polynomial}) = Interval(-Inf, Inf)
mapdomain(::Type{<:Polynomial}, x::AbstractArray) = x

function (p::Polynomial{T})(x::S) where {T,S}
    R = promote_type(T, S)
    length(p) == 0 && return zero(R)
    b = convert(R, p[end])
    @inbounds for i in (lastindex(p) - 1):-1:0
        b = p[i] + x * b
    end
    return b
end

function fromroots(P::Type{<:Polynomial}, r::AbstractVector{T}; var::SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j = 1:n, i = j:-1:1
        c[i + 1] = c[i + 1] - r[j] * c[i]
    end
    return Polynomial(reverse(c), var)
end

function vander(P::Type{<:Polynomial}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    @inbounds for i in 1:n
        A[:, i + 1] = A[:, i] .* x
    end
    return A
end

function integral(p::Polynomial{T}, k::S) where {T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(k)
        return Polynomial([NaN])
    end
    n = length(p)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return Polynomial(a2, p.var)
end

function derivative(p::Polynomial{T}, order::Integer) where {T}
    order < 0 && error("Order of derivative must be non-negative")
    order == 0 && return p
    hasnan(p) && return Polynomial(T[NaN], p.var)
    order > length(p) && return zero(Polynomial{T})

    n = length(p)
    a2 = Vector{T}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return Polynomial(a2, p.var)
end

function companion(p::Polynomial{T}) where T
    d = length(p) - 1
    d < 1 && error("Series must have degree greater than 1")
    d == 1 && return diagm(0 => [-p[0] / p[1]])

    R = eltype(one(T) / p.coeffs[end])
    comp = diagm(-1 => ones(R, d - 1))
    monics = p.coeffs ./ p.coeffs[end]
    comp[:, end] .= -monics[1:d]
    return comp
end

function Base.:+(p1::Polynomial, p2::Polynomial)
    p1.var != p2.var && error("Polynomials must have same variable")
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n]
    return Polynomial(c, p1.var)
end

function Base.:*(p1::Polynomial{T}, p2::Polynomial{S}) where {T,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    n = length(p1) - 1
    m = length(p2) - 1
    R = promote_type(T, S)
    c = zeros(R, m + n + 1)
    for i = 0:n, j = 0:m
        c[i + j + 1] += p1[i] * p2[j]
    end
    return Polynomial(c, p1.var)
end

function Base.divrem(num::Polynomial{T}, den::Polynomial{S}) where {T,S}
    num.var != den.var && error("Polynomials must have same variable")
    n = length(num) - 1
    m = length(den) - 1
    if m == 0 && den[0] ≈ 0 throw(DivideError()) end
    R = typeof(one(T) / one(S))
    P = Polynomial{R}
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

function showterm(io::IO, ::Type{Polynomial{T}}, pj::T, var, j, first::Bool, mimetype) where {T}
    if pj == zero(T) return false end
    pj = printsign(io, pj, first, mimetype)
    if !(pj == one(T) && !(showone(T) || j == 0))   
        printcoefficient(io, pj, j, mimetype)
    end
    printproductsign(io, pj, j, mimetype)
    printexponent(io, var, j, mimetype)
    return true
end
