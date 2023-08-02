# Dense + StandardBasis

"""
    Polynomial{T,X}(coeffs::AbstractVector, [m::Integer = 0], [var = :x])

# Examples:
```jldoctest polynomial
julia> using Polynomials

```
"""
const 𝑃olynomial = MutableDensePolynomial{StandardBasis}
export 𝑃olynomial

_typealias(::Type{P}) where {P<:𝑃olynomial} = "Polynomial"

function ⊗(p::𝑃olynomial{T,X},
           q::𝑃olynomial{S,X}) where {T,S,X}
    # simple convolution
    # This is ⊗(P,p,q) from polynomial standard-basis
    R = promote_type(T,S)
    P = 𝑃olynomial{R,X}

    iszero(p) && return zero(P)
    iszero(q) && return zero(P)

    b₁, b₂ = lastindex(p), lastindex(q)
    b = b₁ + b₂

    z = zero(first(p) * first(q))
    cs = _zeros(p, z, b + 1)

    # convolve and shift order
    @inbounds for (i, pᵢ) ∈ enumerate(p.coeffs)
        for (j, qⱼ) ∈ enumerate(q.coeffs)
            ind = i + j - 1
            cs[ind] = muladd(pᵢ, qⱼ, cs[ind])
        end
    end
    if iszero(last(cs))
        cs = trim_trailing_zeros(cs)
    end
    P(Val(false), cs)
end
