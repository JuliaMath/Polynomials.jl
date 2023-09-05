"""
    SparseVectorPolynomial{T, X}(coeffs::AbstractVector{T}, [var = :x])

```
"""
const SparseVectorPolynomial = MutableSparseVectorPolynomial{StandardBasis}
export SparseVectorPolynomial

_typealias(::Type{P}) where {P<:SparseVectorPolynomial} = "SparseVectorPolynomial"

function evalpoly(x, p::SparseVectorPolynomial)
    tot = zero(p[0]*x)
    for (i, cᵢ) ∈ pairs(p)
        tot = muladd(cᵢ, x^i, tot)
    end
    return tot
end


function Base.:*(p::MutableSparseVectorPolynomial{StandardBasis,T,X},
                 q::MutableSparseVectorPolynomial{StandardBasis,S,X}) where {T,S,X}

    n,m = length(p), length(q)
    R = promote_type(T,S)
    cs = spzeros(T, 1 + n + m)
    @inbounds for (i, aᵢ) ∈ pairs(p)
        for (j, aⱼ) ∈ pairs(q)
            k = i + j
            cs[k + 1] = muladd(aᵢ, aⱼ, cs[k+1])
        end
    end
    MutableSparseVectorPolynomial{StandardBasis,R,X}(cs)
end

function derivative(p:: MutableSparseVectorPolynomial{B,T,X}) where {B<:StandardBasis,T,X}
    n = length(p)
    R = promote_type(T, Int)
    cs = spzeros(R, n)
    @inbounds for (i, aᵢ) ∈ pairs(p)
        cs[i] = i * aᵢ
    end
    MutableSparseVectorPolynomial{StandardBasis,R,X}(cs)
end

function integrate(p:: MutableSparseVectorPolynomial{B,T,X}) where {B<:StandardBasis,T,X}
    n = length(p)
    R = Base.promote_op(/, T, Int)
    cs = spzeros(R, n+1)
    @inbounds for (i, aᵢ) ∈ pairs(p)
        j = i + 1
        iszero(i) && continue
        cs[j + 1] =  aᵢ / j
    end
    MutableSparseVectorPolynomial{StandardBasis,R,X}(cs)
end
