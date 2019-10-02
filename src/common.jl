#=
This file contains all of the abstract definitions for our univariate polynomial type
=#
const SymbolLike = Union{AbstractString,Char,Symbol}

abstract type AbstractPolynomial <: Number end

export AbstractPolynomial,
       roots

using LinearAlgebra: eigvals
import Base: ==, ≈

"""
    roots(::AbstractVector{<:Number}, ::Type{<:AbstractPolynomial}, var=:x)

Construct a polynomial of the given type given the roots.

# Examples
```jldoctest
julia> r = [3, 2]; # (x - 3)(x - 2)

julia> roots(r, Polynomial)
Polynomial(x^2 - 5x + 6)
```
"""
roots(r::AbstractVector{T}, P::Type{<:AbstractPolynomial}, var::SymbolLike = :x) where {T <: Number} = from_roots(P, r, var)

"""
    roots(::AbstractMatrix{<:Number}, ::Type{<:AbstractPolynomial}, var=:x)

Construct a polynomial of the given type using the eigenvalues of the given matrix as the roots.

# Examples
```jldoctest
julia> A = [1 2; 3 4]; # (x - 5.37228)(x + 0.37228)

julia> roots(A, Polynomial)
Polynomial(-1.9999999999999998 - 5.0⋅x + 1.0⋅x^2)
```
"""
roots(A::AbstractMatrix{T}, P::Type{<:AbstractPolynomial}, var::SymbolLike = :x) where {T <: Number} = roots(eigvals(A), P, var)
roots(r, P::Type{<:AbstractPolynomial}, var::SymbolLike = :x) where {T <: Number} = roots(collect(r), P, var)


#=
Inspection
=#
Base.length(p::P) where {P <: AbstractPolynomial} = length(p.coeffs)

#=
Comparisons
=#
==(p1::P, p2::P) where {P <: AbstractPolynomial} = p1.coeffs == p2.coeffs
≈(p1::P, p2::P) where {P <: AbstractPolynomial} = p1.coeffs ≈ p2.coeffs