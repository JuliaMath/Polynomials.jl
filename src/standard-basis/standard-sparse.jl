#const SparsePolynomial = SparseUnivariatePolynomial{StandardBasis} # const is important!
#export SparsePolynomial

function evalpoly(x, p::MutableSparsePolynomial)

    tot = zero(p[0]*x)
    for (i, cᵢ) ∈ p.coeffs
        tot = muladd(cᵢ, x^i, tot)
    end
    return tot
end


function constantterm(p::MutableSparsePolynomial{B,T,X}) where {B,T,X}
    get(p.coeffs, 0, zero(T))
end

function ⊗(p::MutableSparsePolynomial{StandardBasis,T,X},
           q::MutableSparsePolynomial{StandardBasis,S,X}) where {T,S,X}
    # simple convolution
    R = promote_type(T,S)
    P = MutableSparsePolynomial{StandardBasis,R,X}

    z = zero(R)
    cs = Dict{Int, R}()

    @inbounds for (i, pᵢ) ∈ pairs(p)
        for (j, qⱼ) ∈ pairs(q)
            cᵢⱼ = get(cs, i+j, z)
            cs[i+j] = muladd(pᵢ, qⱼ, cᵢⱼ)
        end
    end

    P(cs)
end
