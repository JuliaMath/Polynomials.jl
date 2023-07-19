# Dense + StandardBasis
#const Polynomial = MutableDensePolynomial{StandardBasis} # const is important!
#export Polynomial

#LaurentPolynomial = MutableDensePolynomial{StandardBasis}
#export LaurentPolynomial

function evalpoly(c, p::MutableDensePolynomial{StandardBasis,T,X}) where {T,X}
    iszero(p) && return zero(T)*zero(c)
    EvalPoly.evalpoly(c, p.coeffs) * c^p.order
end


# function isconstant(p:: MutableDensePolynomial{StandardBasis})
#     firstindex(p) != 0 && return false
#     i = findlast(!iszero, p.coeffs)
#     i == nothing && return true
#     i == 1 && return true
# end

# scalar add
function scalar_add(c::S, p:: MutableDensePolynomial{StandardBasis,T,X}) where {S, T, X}
    R = promote_type(T,S)
    P =  MutableDensePolynomial{StandardBasis,R,X}

    iszero(p) && return P([c], 0)
    iszero(c) && return convert(P, p)

    a,b = firstindex(p), lastindex(p)
    a′ = min(0,a)
    cs = _zeros(p, zero(first(p.coeffs)+c), length(a′:b))
    o = offset(p) + a - a′
    for (i, cᵢ) ∈ pairs(p)
        cs[i+o] = cᵢ
    end
    cs[0+o] += c
    iszero(last(cs)) && (cs = trim_trailing_zeros(cs))
    P(Val(false), cs, a′)
end



function ⊗(p:: MutableDensePolynomial{StandardBasis,T,X},
           q:: MutableDensePolynomial{StandardBasis,S,X}) where {T,S,X}
    # simple convolution
    R = promote_type(T,S)
    P =  MutableDensePolynomial{StandardBasis,R,X}

    iszero(p) && return zero(P)
    iszero(q) && return zero(P)

    a₁, a₂ = firstindex(p), firstindex(q)
    b₁, b₂ = lastindex(p), lastindex(q)
    a, b = a₁ + a₂, b₁ + b₂

    z = zero(first(p) * first(q))
    cs = _zeros(p, z, length(a:b))

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
    P(Val(false), cs, a)
end

# function ⊗(p:: MutableDensePolynomial{StandardBasis}{T,X},
#            q:: MutableDensePolynomial{StandardBasis}{S,X}) where {T,S,X}
#     # simple convolution
#     R = promote_type(T,S)
#     P =  MutableDensePolynomial{StandardBasis}{R,X}

#     iszero(p) && return zero(P)
#     iszero(q) && return zero(P)

#     a₁, a₂ = firstindex(p), firstindex(q)
#     b₁, b₂ = lastindex(p), lastindex(q)
#     a, b = a₁ + a₂, b₁ + b₂

#     z = zero(first(p) * first(q))
#     cs = _zeros(p, z, length(a:b))

#     # convolve and shift order
#     @inbounds for (i, pᵢ) ∈ enumerate(p.coeffs)
#         for (j, qⱼ) ∈ enumerate(q.coeffs)
#             ind = i + j - 1
#             cs[ind] += pᵢ * qⱼ
#         end
#     end

#     iszero(last(cs)) && chop_right!(cs)
#     P(Val(false), cs, a)
# end
