export Polynomial

"""
    Polynomial{T, X}(coeffs::AbstractVector{T}, [var = :x])

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
julia> using Polynomials

julia> Polynomial([1, 0, 3, 4])
Polynomial(1 + 3*x^2 + 4*x^3)

julia> Polynomial([1, 2, 3], :s)
Polynomial(1 + 2*s + 3*s^2)

julia> one(Polynomial)
Polynomial(1.0)
```
"""
struct Polynomial{T, X} <: StandardBasisPolynomial{T, X}
    coeffs::Vector{T}
    function Polynomial{T, X}(coeffs::AbstractVector{S}) where {T, X, S}
        if Base.has_offset_axes(coeffs)
            @warn "ignoring the axis offset of the coefficient vector"
        end
        N = findlast(!iszero, coeffs)
        isnothing(N) && return new{T,X}(T[])
        cs = T[coeffs[i] for i ∈ firstindex(coeffs):N]
        new{T,X}(cs)
    end
    # non-copying alternative assuming !iszero(coeffs[end])
    function Polynomial{T, X}(checked::Val{false}, coeffs::AbstractVector{T}) where {T, X}
        new{T, X}(coeffs)
    end
end

@register Polynomial

# some speedups using non-copying constructor
Base.zero(::Type{P}) where {T, X, P<:Polynomial{T,X}} =
    Polynomial{T,X}(Val(false), T[])

# scalar +,* faster  than standard-basis/common versions as it avoids a copy
function Base.:+(p::P, c::S) where {T, X, P <: Polynomial{T, X}, S<:Number}
    if isempty(p.coeffs)
        as = promote_type(S,T)[c]
    else
        R = typeof(c + p.coeffs[1])
        as = R[isone(i) ? cᵢ + c : cᵢ for (i,cᵢ) ∈ pairs(p.coeffs)]
    end
    _polynomial(p, as)
end

# scalar * is a bit faster (2 to 3 allocations), isave allocations on copying
function scalar_mult(p::P, c::S) where {T, X, P <: Polynomial{T,X} , S <: Number}
    return _polynomial(p, coeffs(p) .* (c,))
end

function scalar_mult(c::S, p::P) where {T, X, P <: Polynomial{T,X} , S <: Number}
    _polynomial(p,  (c,).* coeffs(p))
end
# much faster than Polynomial(as, X)
function _polynomial(p::P, as)  where {T, X, P <: Polynomial{T,X}}
    R = eltype(as)
    Q = Polynomial{R, X}
    isempty(as) && return Q(Val(false), as)
    iszero(as[end]) ? Q(as) : Q(Val(false), as)
end


# implement, as not copying speeds up multiplication by a factor of 2 or so
# over the default
function Base.:+(p1::P1, p2::P2) where {T,X, P1<:Polynomial{T,X},
                                        S,   P2<:Polynomial{S,X}}
    n1, n2 = length(p1), length(p2)
    R = promote_type(T,S)
    Q = Polynomial{R,X}
    if n1 == n2
        cs = ⊕(P1, p1.coeffs, p2.coeffs)
        isempty(cs) && return Q(Val(false), cs)
        return iszero(cs[end]) ? Q(cs) : Q(Val(false), cs)
    elseif n1 > n2
        cs = ⊕(P1, p1.coeffs, p2.coeffs)
    else
        cs = ⊕(P1, p2.coeffs, p1.coeffs)
    end

    Q(Val(false), cs)
end

# redundant, a bit faster
function Base.:*(p::P, q::P) where {T <: Number,X, P<:Polynomial{T,X}}
    c = fastconv(p.coeffs, q.coeffs)
    (isempty(c) || !iszero(last(c))) && return P(Val(false), c)
    return P(c)
end
