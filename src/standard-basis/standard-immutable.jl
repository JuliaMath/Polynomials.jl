#const ImmutablePolynomial = ImmutableDensePolynomial{StandardBasis}
#export ImmutablePolynomial

function evalpoly(x, p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,N}
    N == 0 && return zero(T) * zero(x)
    z = zero(x * zero(p[0]))
    typeof(z)(EvalPoly.evalpoly(x, p.coeffs))
end

# Padded vector sum of two tuples assuming N ≥ M
@generated function tuple_sum(p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}

    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(p1[$i] + p2[$i])
    end
    for i in M+1:N
        exprs[i] =:(p1[$i])
    end

    return quote
        Base.@_inline_meta
        #Base.@inline
        tuple($(exprs...))
    end

end



function scalar_add(p::ImmutableDensePolynomial{B,T,X,N}, c::S) where {B<:StandardBasis,T,X,S,N}
    R = promote_type(T,S)
    P = ImmutableDensePolynomial{B,R,X}
    iszero(c) && return P{N}(convert(NTuple{N,R}, p.coeffs))
    N == 0 && return P{1}(NTuple{1,R}(c))
    N == 1 && return P{N}((p[0]+c,))

    #cs = tuple_sum(convert(NTuple{N,R}, p.coeffs), NTuple{1,R}(c))
    cs = _tuple_combine(+, convert(NTuple{N,R}, p.coeffs), NTuple{1,R}(c))
    q = P{N}(cs)

    return q


    P = ImmutableDensePolynomial{B,R,X}
    iszero(N) && return P{1}((c,))

    xs = convert(NTuple{N,R}, p.coeffs)
    @set! xs[1] = xs[1] + c
    P{N}(xs)
end


function XXscalar_add(p::ImmutableDensePolynomial{B,T,X,N}, c::S) where {B,T,X,S,N}
    R = promote_type(T,S)
    P = ImmutableDensePolynomial{B,R,X}
    iszero(N) && return P{1}((c,))

    xs = convert(NTuple{N,R}, p.coeffs)
    @set! xs[1] = xs[1] + c
    P{N}(xs)
end

# return N*M
function ⊗(p::ImmutableDensePolynomial{StandardBasis,T,X,N},
           q::ImmutableDensePolynomial{StandardBasis,S,X,M}) where {T,S,X,N,M}

    # simple convolution
    R = promote_type(T,S)
    P = ImmutableDensePolynomial{StandardBasis,R,X}

    (iszero(N) || iszero(M)) && return zero(P)

    cs = fastconv(p.coeffs, q.coeffs)
    P{N+M-1}(cs)
end

## Static size of product makes generated functions  a good choice
## from https://github.com/tkoolen/StaticUnivariatePolynomials.jl/blob/master/src/monomial_basis.jl
## convolution of two tuples
@generated function fastconv(p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}
    P = M + N - 1
    exprs = Any[nothing for i = 1 : P]
    for i in 1 : N
        for j in 1 : M
            k = i + j - 1
            if isnothing(exprs[k])
                exprs[k] = :(p1[$i] * p2[$j])
            else
                exprs[k] = :(muladd(p1[$i], p2[$j], $(exprs[k])))
            end
        end
    end

    return quote
        Base.@_inline_meta # 1.8 deprecation
        tuple($(exprs...))
    end

end


function differentiate(p::ImmutableDensePolynomial{StandardBasis,T,X,N}) where {T,X,N}
    N == 0 && return 1p
    cs = ntuple(i -> i*p.coeffs[i+1], Val(N-1))
    R = eltype(cs)
    ImmutableDensePolynomial{StandardBasis,R,X,N-1}(cs)
end


function integrate(p::ImmutableDensePolynomial{StandardBasis,T,X,N}) where {T,X,N}
    cs = ntuple(i -> i > 1 ? p.coeffs[i-1]/(i-1) : zero(T)/1, Val(N+1))
    R = eltype(cs)
    ImmutableDensePolynomial{StandardBasis,R,X,N+1}(cs)
end
