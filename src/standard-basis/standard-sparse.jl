## Standard basis + sparse storage

"""
    SparsePolynomial{T, X}(coeffs::Dict{Int,T})

Polynomials in the standard basis backed by a dictionary holding the
non-zero coefficients. For polynomials of high degree, this might be
advantageous.

# Examples:

```jldoctest
julia> using Polynomials

julia> P  = SparsePolynomial;

julia> p, q = P([1,2,3]), P([4,3,2,1])
(SparsePolynomial(1 + 2*x + 3*x^2), SparsePolynomial(4 + 3*x + 2*x^2 + x^3))

julia> p + q
SparsePolynomial(5 + 5*x + 5*x^2 + x^3)

julia> p * q
SparsePolynomial(4 + 11*x + 20*x^2 + 14*x^3 + 8*x^4 + 3*x^5)

julia> p + 1
SparsePolynomial(2 + 2*x + 3*x^2)

julia> q * 2
SparsePolynomial(8 + 6*x + 4*x^2 + 2*x^3)

julia> p = Polynomials.basis(P, 10^9) - Polynomials.basis(P,0) # also P(Dict(0=>-1, 10^9=>1))
SparsePolynomial(-1.0 + 1.0*x^1000000000)

julia> p(1)
0.0
```

!!! note
    `SparsePolynomial` is an alias for `MutableSparsePolynomial{StandardBasis}`.

"""
const SparsePolynomial = MutableSparsePolynomial{StandardBasis} # const is important!
export SparsePolynomial

_typealias(::Type{P}) where {P<:SparsePolynomial} = "SparsePolynomial"

function evalpoly(x, p::MutableSparsePolynomial)

    tot = zero(p[0]*x)
    for (i, cᵢ) ∈ p.coeffs
        tot = muladd(cᵢ, x^i, tot)
    end
    return tot
end

# much faster than default
function scalar_add(c::S, p::MutableSparsePolynomial{B,T,X}) where {B<:StandardBasis,T,X,S}
    c₀ = c + p[0]
    R = eltype(c₀)
    P = MutableSparsePolynomial{B,R,X}
    D = convert(Dict{Int, R}, copy(p.coeffs))
    if iszero(c₀)
        delete!(D,0)
    else
        @inbounds D[0] = c₀
    end
    return P(Val(false), D)
end


function ⊗(p::MutableSparsePolynomial{StandardBasis,T,X},
           q::MutableSparsePolynomial{StandardBasis,S,X}) where {T,S,X}

    # simple convolution same as ⊗(Polynomial, p.coeffs, q.coeffs) (DRY)
    R = promote_type(T,S)
    P = MutableSparsePolynomial{StandardBasis,R,X}

    z = zero(R)
    cs = Dict{Int, R}()

    for (i, pᵢ) ∈ pairs(p)
        for (j, qⱼ) ∈ pairs(q)
            cᵢⱼ = get(cs, i+j, z)
            val = muladd(pᵢ, qⱼ, cᵢⱼ)
            iszero(val) && continue
            @inbounds cs[i+j] = val
        end
    end
    P(Val(false), cs)
end


# sparse
function derivative(p:: MutableSparsePolynomial{B,T,X}) where {B<:StandardBasis,T,X}
    N = lastindex(p) - firstindex(p) + 1
    R = promote_type(T, Int)
    P = ⟒(p){R,X}
    hasnan(p) && return  P(zero(T)/zero(T)) # NaN{T}
    iszero(p) && return zero(P)

    d = Dict{Int,R}()
    for (i, pᵢ) ∈ pairs(p)
        iszero(i) && continue
        d[i-1] = i*pᵢ
    end
    return P(d)
end

function integrate(p:: MutableSparsePolynomial{B,T,X}) where {B<:StandardBasis,T,X}

    R = Base.promote_op(/, T, Int)
    P = MutableSparsePolynomial{B,R,X}
    hasnan(p) && return ⟒(P)(NaN)
    iszero(p) && return zero(p)/1

    d = Dict{Int, R}()
    for (i, pᵢ) ∈ pairs(p.coeffs)
        i == -1 && throw(ArgumentError("Can't integrate Laurent polynomial with  `x⁻¹` term"))
        cᵢ₊₁ = pᵢ/(i+1)
        !iszero(cᵢ₊₁) && (d[i+1] = cᵢ₊₁)
    end
    return P(d)
end
