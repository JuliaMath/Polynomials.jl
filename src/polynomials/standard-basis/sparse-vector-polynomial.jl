"""
    SparseVectorPolynomial{T, [X]}(coeffs::AbstractVector{T}, [var = :x])

A polynomial using a sparse vector for holding the coefficients.

Expects `zero(T)` to be defined.

See also [SparsePolynomial](@ref) which uses a dictionary to store the coefficients. Compared to that, some basic operations are faster, some not so (addition).
```
"""
const SparseVectorPolynomial = MutableSparseVectorPolynomial{StandardBasis}
#export SparseVectorPolynomial

_typealias(::Type{P}) where {P<:SparseVectorPolynomial} = "SparseVectorPolynomial"

# generally faster to use zip(findnz(p.coeffs)...) with an offset than pairs(p), which sorts
function evalpoly(x, p::SparseVectorPolynomial)
    tot = zero(p[0]*x)
    for (i, cᵢ) ∈ zip(findnz(p.coeffs)...)
        j = i - 1
        tot = muladd(cᵢ, x^j, tot)
    end
    return tot
end

function scalar_add(c::S, p::P) where {B<:StandardBasis, S, T, X, P<:MutableSparseVectorPolynomial{B,T,X}}
    R = promote_type(S,T)
    Q = MutableSparseVectorPolynomial{B,R,X}
    isempty(p.coeffs) && return Q([c])

    cs = similar(p.coeffs, R)
    copy!(cs, p.coeffs)
    cs[1] += c
    Q(cs)

end

function Base.:*(p::MutableSparseVectorPolynomial{StandardBasis,T,X},
                 q::MutableSparseVectorPolynomial{StandardBasis,S,X}) where {T,S,X}

    n,m = length(p), length(q)
    R = promote_type(T,S)
    cs = spzeros(R, 1 + n + m)

    @inbounds for (i, aᵢ) ∈ zip(findnz(p.coeffs)...) #pairs(p)
        for (j, aⱼ) ∈ zip(findnz(q.coeffs)...) #pairs(q)
            k = i + j - 2
            cs[k + 1] = muladd(aᵢ, aⱼ, cs[k+1])
        end
    end

    MutableSparseVectorPolynomial{StandardBasis,R,X}(cs)
end

function derivative(p:: MutableSparseVectorPolynomial{B,T,X}) where {B<:StandardBasis,T,X}
    n = length(p)
    R = promote_type(T, Int)
    P = MutableSparseVectorPolynomial{StandardBasis,R,X}
    iszero(n) && return zero(P)

    cs = spzeros(R, n)

    @inbounds for (i, aᵢ) ∈ zip(findnz(p.coeffs)...) #pairs(p)
        !isfinite(aᵢ) && isnan(aᵢ) && return P([NaN])
        j = i - 1
        iszero(j) && continue
        cs[j] = j * aᵢ
    end

    P(cs)
end

function integrate(p::MutableSparseVectorPolynomial{B,T,X}) where {B<:StandardBasis,T,X}
    n = length(p)
    R = Base.promote_op(/, T, Int)
    P = MutableSparseVectorPolynomial{StandardBasis,R,X}
    iszero(n) && return zero(P)
    cs = spzeros(R, n+1)

    @inbounds for (i, aᵢ) ∈ zip(findnz(p.coeffs)...) #pairs(p)
        !isfinite(aᵢ) && isnan(aᵢ) && return P([NaN])
        j = i# + 1
        cs[j + 1] =  aᵢ / j
    end

    P(cs)
end
