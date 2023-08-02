"""
    PnPolynomial{T,X}(coeffs::Vector{T})

Construct a polynomial in `P_n` (or `Πₙ`), the collection of polynomials in the
standard basis of degree `n` *or less*, using a vector of length
`N+1`.

* Unlike other polynomial types, this type allows trailing zeros in the coefficient vector. Call `chop!` to trim trailing zeros if desired.
* Unlike other polynomial types, this does not copy the coefficients on construction
* Unlike other polynomial types, this type broadcasts like a vector for in-place vector operations (scalar multiplication, polynomial addition/subtraction of the same size)
* The inplace method `mul!(pq, p, q)` is defined to use precallocated storage for the product of `p` and `q`

This type is useful for reducing copies and allocations in some algorithms.

"""
const PnPolynomial = MutableDenseViewPolynomial{StandardBasis}
_typealias(::Type{P}) where {P<:PnPolynomial} = "PnPolynomial"

function ⊗(p:: PnPolynomial{T,X},
           q:: PnPolynomial{S,X}) where {T,S,X}
    R = promote_type(T,S)
    P =  PnPolynomial{R,X}

    N,M = length(p), length(q)
    iszero(N) && return zero(P)
    iszero(M) && return zero(P)

    cs = Vector{R}(undef, N+M-1)
    pq = P(cs)
    mul!(pq, p, q)
    return pq
end

# This is for standard basis
function LinearAlgebra.mul!(pq::PnPolynomial, p::PnPolynomial, q::PnPolynomial)
    m,n = length(p)-1, length(q)-1
    cs = pq.coeffs
    cs .= 0
    @inbounds for (i, pᵢ) ∈ enumerate(p.coeffs)
        for (j, qⱼ) ∈ enumerate(q.coeffs)
            ind = i + j - 1
            cs[ind] = muladd(pᵢ, qⱼ, cs[ind])
        end
    end
    nothing
end
