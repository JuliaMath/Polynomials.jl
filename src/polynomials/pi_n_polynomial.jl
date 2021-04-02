""" 

A polynomial in Πₘ, meaning we keep n+1 coefficient vectors, though possibly the tail ones are zeros.
"""


"""
    ΠₙPolynomial{T,X}(coeffs::Vector{T})

Construct a polynomial in `Πₙ`, the collection of polynomials of degree `n` or less using a vector of length `N+1`.

* Unlike other polynomial types, this does not copy the coefficients on construction
* Unlike other polynomial types, this type broadcasts like a vector for in-place vector operations (scalare multiplication, polynomial addition/subtraction of the same size)

"""
struct ΠₙPolynomial{T,X} <: Polynomials.StandardBasisPolynomial{T, X}
    coeffs::Vector{T}
    function ΠₙPolynomial{T, X}(coeffs::AbstractVector{T}) where {T, X}
        N = length(coeffs) - 1
        new{T,X}(coeffs) # NO CHECK on trailing zeros
    end
end



Polynomials.@register ΠₙPolynomial

# change broadcast semantics
Base.broadcastable(p::ΠₙPolynomial) = p.coeffs;
Base.ndims(::Type{<:ΠₙPolynomial}) = 1
Base.copyto!(p::ΠₙPolynomial, x) = copyto!(p.coeffs, x);

maxdegree(p::ΠₙPolynomial{T,X}) where {T,X} = length(p)-1
function Polynomials.degree(p::ΠₙPolynomial)
    i = findlast(!iszero, p.coeffs)
    i == nothing && return -1
    i - 1
end

# pre-allocated multiplication
function LinearAlgebra.mul!(pq, p::ΠₙPolynomial{T,X}, q) where {T,X}
    m,n = maxdegree(p), maxdegree(q)    
    pq.coeffs .= zero(T)
    for i ∈ 0:m
        for j ∈ 0:n
            k = i + j
            pq.coeffs[1+k] += p.coeffs[1+i] * q.coeffs[1+j]
        end
    end
    nothing
end

# also just p .*= -1
function Base.:-(p::ΠₙPolynomial{T,X}) where {T,X}
    for i ∈ eachindex(p.coeffs)
        p.coeffs[i] *= -1
    end
    p
end



