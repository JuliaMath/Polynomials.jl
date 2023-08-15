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

# used by multroot
function LinearAlgebra.mul!(pq::PnPolynomial, p::PnPolynomial{T,X}, q::PnPolynomial) where {T,X}
    m,n = length(p)-1, length(q)-1
    pq.coeffs .= zero(T)
    for i ∈ 0:m
        for j ∈ 0:n
            k = i + j
            @inbounds pq.coeffs[1+k] = muladd(p.coeffs[1+i], q.coeffs[1+j], pq.coeffs[1+k])
        end
    end
    nothing
end
