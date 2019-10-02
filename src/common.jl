#=
This file contains all of the abstract definitions for our univariate polynomial type
=#
const SymbolLike = Union{AbstractString,Char,Symbol}

abstract type AbstractPolynomial <: Number end

export AbstractPolynomial,
       fromroots

using LinearAlgebra: eigvals
import Base: ==, ≈

"""
    fromroots(::AbstractVector{<:Number}, var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractVector{<:Number}, var=:x)

Construct a polynomial of the given type given the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest
julia> r = [3, 2]; # (x - 3)(x - 2)

julia> fromroots(r)
Polynomial(x^2 - 5x + 6)
```
"""
fromroots(P::Type{<:AbstractPolynomial}, r::AbstractVector{T}, var::SymbolLike = :x) where {T <: Number} = _fromroots(P, r, var)
fromroots(r::AbstractVector{<:Number}, var::SymbolLike = :x) = fromroots(Polynomial, r, var)
fromroots(r, var::SymbolLike = :x) = fromroots(collect(r), var)

"""
    fromroots(::AbstractMatrix{<:Number}, var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractMatrix{<:Number}, var=:x)

Construct a polynomial of the given type using the eigenvalues of the given matrix as the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest
julia> A = [1 2; 3 4]; # (x - 5.37228)(x + 0.37228)

julia> fromroots(A)
Polynomial(-1.9999999999999998 - 5.0⋅x + 1.0⋅x^2)
```
"""
fromroots(P::Type{<:AbstractPolynomial}, A::AbstractMatrix{T}, var::SymbolLike = :x) where {T <: Number} = fromroots(P, eigvals(A), var)
fromroots(A::AbstractMatrix{T}, var::SymbolLike = :x) where {T <: Number} = fromroots(Polynomial, eigvals(A), var)

fit(x, y, P::Type{<:AbstractPolynomial}) = fit(collect(x), collect(y), P)
function fit(x::AbstractVector, y::AbstractVector, P::Type{<:AbstractPolynomial})
    x = scale_to_domain(P, x)
    vand = vander(P, x)
    coeffs = pinv(vand)
    return P(coeffs)
end

#=
Inspection
=#
Base.length(p::P) where {P <: AbstractPolynomial} = length(p.coeffs)

#=
Comparisons
=#
==(p1::P, p2::P) where {P <: AbstractPolynomial} = p1.coeffs == p2.coeffs
≈(p1::P, p2::P) where {P <: AbstractPolynomial} = p1.coeffs ≈ p2.coeffs