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
const 𝑃olynomial = MutableDensePolynomial{StandardBasis}
const Polynomial = 𝑃olynomial
export Polynomial

_typealias(::Type{P}) where {P<:𝑃olynomial} = "𝑃olynomial"
