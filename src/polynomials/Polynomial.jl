export Polynomial

"""
    Polynomial{T<:Number, X}(coeffs::AbstractVector{T}, [var = :x])

Construct a polynomial from its coefficients `coeffs`, lowest order first, optionally in
terms of the given variable `var` which may be a character, symbol, or a string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct this through
`Polynomial([a_0, a_1, ..., a_n])`.

The usual arithmetic operators are overloaded to work with polynomials as well as
with combinations of polynomials and scalars. However, operations involving two
polynomials of different variables causes an error except those involving a constant polynomial.

!!! note
    `Polynomial` is not axis-aware, and it treats `coeffs` simply as a list of coefficients with the first 
    index always corresponding to the constant term. In order to use the axis of `coeffs` as exponents, 
    consider using a [`LaurentPolynomial`](@ref) or possibly a [`SparsePolynomial`](@ref).

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
struct Polynomial{T <: Number, X} <: StandardBasisPolynomial{T, X}
    coeffs::Vector{T}
    function Polynomial{T, X}(coeffs::AbstractVector{S}) where {T <: Number, X, S}
        if Base.has_offset_axes(coeffs)
            @warn "ignoring the axis offset of the coefficient vector"
        end
        N = findlast(!iszero, coeffs)
        isnothing(N) && return new{T,X}(zeros(T,1))
        cs = T[coeffs[i] for i ∈ firstindex(coeffs):N]
        new{T,X}(cs)
    end
    # non-copying alternative assuming !iszero(coeffs[end])
    function Polynomial{T, X}(checked::Val{false}, coeffs::AbstractVector{T}) where {T <: Number, X}
        new{T, X}(coeffs)
    end
end

@register Polynomial

"""
    (p::Polynomial)(x)

Evaluate the polynomial using [Horner's Method](https://en.wikipedia.org/wiki/Horner%27s_method), also known as synthetic division, as implemented in `evalpoly` of base `Julia`.

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


# scalar _,* faster  than standard-basis/common versions
function Base.:+(p::P, c::S) where {T, X, P <: Polynomial{T, X}, S<:Number}
    R = promote_type(T, S)
    Q = Polynomial{R,X}
    as = convert(Vector{R}, copy(coeffs(p)))
    as[1] += c
    iszero(as[end]) ? Q(as) : Q(Val(false), as)
end

function Base.:*(p::P, c::S) where {T, X, P <: Polynomial{T,X} , S <: Number}
    R = promote_type(T,S)
    Q = Polynomial{R, X}
    as = R[aᵢ * c for aᵢ ∈ coeffs(p)]
    iszero(as[end]) ? Q(as) : Q(Val(false), as)
end

function Base.:+(p1::Polynomial{T}, p2::Polynomial{S}) where {T, S}
    isconstant(p1) && return p2 + p1[0]
    isconstant(p2) && return p1 + p2[0]
    X, Y = indeterminate(p1), indeterminate(p2)
    X != Y && error("Polynomials must have same variable")
    n1, n2 = length(p1), length(p2)
    R = promote_type(T,S)

    c = zeros(R, max(n1, n2))
    if n1 >= n2
        c .= p1.coeffs
        for i in eachindex(p2.coeffs)
            c[i] += p2.coeffs[i]
        end
    else
        c .= p2.coeffs
        for i in eachindex(p1.coeffs)
            c[i] += p1.coeffs[i]
            end
    end

    Q = Polynomial{R,X}
    return iszero(c[end]) ? Q(c) : Q(Val(false), c)

end


function Base.:*(p1::Polynomial{T}, p2::Polynomial{S}) where {T,S}

    n, m = length(p1)-1, length(p2)-1 # not degree, so pNULL works
    X, Y = indeterminate(p1), indeterminate(p2)
    R = promote_type(T, S)
    if n > 0 && m > 0
        X != Y && error("Polynomials must have same variable")
        c = zeros(R, m + n + 1)
        for i in 0:n, j in 0:m
            @inbounds c[i + j + 1] += p1[i] * p2[j]
        end
        Q = Polynomial{R,X}
        return iszero(c[end]) ? Q(c) : Q(Val(false), c)
    elseif n <= 0
        return Polynomial{R, Y}(p2.coeffs * p1[0])
    else
        return Polynomial{R, X}(p1.coeffs * p2[0])
    end

end
