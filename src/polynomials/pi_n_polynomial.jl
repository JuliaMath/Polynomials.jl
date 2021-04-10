"""
    PnPolynomial{T,X}(coeffs::Vector{T})

Construct a polynomial in `P_n` (or `Πₙ`), the collection of polynomials in the
standard basis of degree `n` *or less*, using a vector of length
`N+1`.

* Unlike other polynomial types, this type allows trailing zeros in the coefficient vector. Call `chop!` to trim trailing zeros if desired.
* Unlike other polynomial types, this does not copy the coefficients on construction
* Unlike other polynomial types, this type broadcasts like a vector for in-place vector operations (scalare multiplication, polynomial addition/subtraction of the same size)
* The method inplace `mul!(pq, p, q)` is defined to use precallocated storage for the product of `p` and `q`

This type is useful for reducing copies and allocations in some algorithms.

"""
struct PnPolynomial{T,X} <: Polynomials.StandardBasisPolynomial{T, X}
    coeffs::Vector{T}
    function PnPolynomial{T, X}(coeffs::AbstractVector{T}) where {T, X}
        N = length(coeffs) - 1
        new{T,X}(coeffs) # NO CHECK on trailing zeros
    end
end



Polynomials.@register PnPolynomial

# change broadcast semantics
Base.broadcastable(p::PnPolynomial) = p.coeffs;
Base.ndims(::Type{<:PnPolynomial}) = 1
Base.copyto!(p::PnPolynomial, x) = copyto!(p.coeffs, x);

function Polynomials.degree(p::PnPolynomial)
    i = findlast(!iszero, p.coeffs)
    i == nothing && return -1
    i - 1
end

# pre-allocated multiplication
function LinearAlgebra.mul!(pq, p::PnPolynomial{T,X}, q) where {T,X}
    m,n = length(p)-1, length(q)-1
    pq.coeffs .= zero(T)
    for i ∈ 0:m
        for j ∈ 0:n
            k = i + j
            pq.coeffs[1+k] += p.coeffs[1+i] * q.coeffs[1+j]
        end
    end
    nothing
end

