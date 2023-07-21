function evalpoly(x, p::MutableSparsePolynomial)

    tot = zero(p[0]*x)
    for (i, cᵢ) ∈ p.coeffs
        tot = muladd(cᵢ, x^i, tot)
    end
    return tot
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

## ---
const SparsePolynomial = MutableSparsePolynomial{StandardBasis} # const is important!
export SparsePolynomial
