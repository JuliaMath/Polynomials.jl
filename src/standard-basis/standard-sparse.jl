#const SparsePolynomial = SparseUnivariatePolynomial{StandardBasis} # const is important!
#export SparsePolynomial

function evalpoly(x, p::SparseUnivariatePolynomial)

    tot = zero(p[0]*x)
    for (i, cᵢ) ∈ p.coeffs
        tot = muladd(cᵢ, x^i, tot)
    end
    return tot
end


function constantterm(p::SparseUnivariatePolynomial{B,T,X}) where {B,T,X}
    get(p.coeffs, 0, zero(T))
end

function ⊗(p::SparseUnivariatePolynomial{StandardBasis,T,X},
           q::SparseUnivariatePolynomial{StandardBasis,S,X}) where {T,S,X}
    # simple convolution
    R = promote_type(T,S)
    P = SparseUnivariatePolynomial{StandardBasis,R,X}

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
