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
```@meta
DocTestSetup = quote
    using Polynomials
end
```

```jldoctest
julia> Polynomial([1, 0, 3, 4])
Polynomial(1 + 3*x^2 + 4*x^3)

julia> Polynomial([1, 2, 3], :s)
Polynomial(1 + 2*s + 3*s^2)

julia> one(Polynomial)
Polynomial(1.0)
```
"""
struct Polynomial{T <: Number} <: StandardBasisPolynomial{T}
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


"""
    (p::Polynomial)(x)

Evaluate the polynomial using [Horner's Method](https://en.wikipedia.org/wiki/Horner%27s_method), also known as synthetic division.

```@meta
DocTestSetup = quote
    using Polynomials
end
```

# Examples
```jldoctest
julia> p = Polynomial([1, 0, 3])
Polynomial(1 + 3*x^2)

julia> p(0)
1

julia> p.(0:3)
4-element Array{Int64,1}:
  1
  4
 13
 28
```
"""
(p::Polynomial{T})(x::S) where {T,S} = evalpoly(x, coeffs(p))



function Base.:+(p1::Polynomial{T}, p2::Polynomial{S}) where {T, S}
    R = promote_type(T,S)
    p1.var != p2.var && error("Polynomials must have same variable")

    n1, n2 = length(p1), length(p2)
    c = [p1[i] + p2[i] for i = 0:max(n1, n2)]
    return Polynomial(c, p1.var)
end


function Base.:+(p::Polynomial{T}, c::S) where {T,S<:Number}
    U = promote_type(T, S)
    q = copy(p)
    p2 = U == S ? q : convert(Polynomial{U}, q)
    p2[0] += c
    return p2
end

function Base.:*(p1::Polynomial{T}, p2::Polynomial{S}) where {T,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    n,m = length(p1)-1, length(p2)-1 # not degree, so pNULL works
    R = promote_type(T, S)
    c = zeros(R, m + n + 1)
    for i in 0:n, j in 0:m
        @inbounds c[i + j + 1] += p1[i] * p2[j]
    end
    return Polynomial(c, p1.var)
end
