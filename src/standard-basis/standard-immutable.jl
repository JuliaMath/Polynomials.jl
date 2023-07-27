## Immutable dense / standard basis specific polynomial code
ImmutablePolynomial = ImmutableDensePolynomial{StandardBasis}
export ImmutablePolynomial


evalpoly(x, p::ImmutableDensePolynomial{B,T,X,0}) where {B<:StandardBasis,T,X} = zero(T)*zero(x)
evalpoly(x, p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,N} = EvalPoly.evalpoly(x, p.coeffs)

constantterm(p::ImmutableDensePolynomial{B,T,X,0}) where {B <: StandardBasis,T,X} = zero(T)
constantterm(p::ImmutableDensePolynomial{B,T,X,N}) where {B <: StandardBasis,T,X,N} = p.coeffs[1]

scalar_add(p::ImmutableDensePolynomial{B,T,X,0}, c::S) where {B<:StandardBasis,T,X,S} =
    ImmutableDensePolynomial{B,promote_type(T,S),X,1}((c,))

function scalar_add(p::ImmutableDensePolynomial{B,T,X,N}, c::S) where {B<:StandardBasis,T,X,S,N}
    R = promote_type(T,S)
    P = ImmutableDensePolynomial{B,R,X}
    iszero(c) && return P{N}(convert(NTuple{N,R}, p.coeffs))
    N == 0 && return P{1}(NTuple{1,R}(c))
    N == 1 && return P{N}((p[0]+c,))

    cs = _tuple_combine(+, convert(NTuple{N,R}, p.coeffs), NTuple{1,R}((c,)))
    q = P{N}(cs)

    return q


    P = ImmutableDensePolynomial{B,R,X}
    iszero(N) && return P{1}((c,))

    xs = convert(NTuple{N,R}, p.coeffs)
    @set! xs[1] = xs[1] + c
    P{N}(xs)
end


# return N*M
# intercept promotion call
function Base.:*(p::ImmutableDensePolynomial{StandardBasis,T,X,N},
                 q::ImmutableDensePolynomial{StandardBasis,S,X,M}) where {T,S,X,N,M}
    ⊗(p,q)
end

⊗(p::ImmutableDensePolynomial{B,T,X,0},
  q::ImmutableDensePolynomial{B,S,X,M}) where {B<:StandardBasis,T,S,X,M} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})
⊗(p::ImmutableDensePolynomial{B,T,X,N},
  q::ImmutableDensePolynomial{B,S,X,0}) where {B<:StandardBasis,T,S,X,N} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})
⊗(p::ImmutableDensePolynomial{B,T,X,0},
  q::ImmutableDensePolynomial{B,S,X,0}) where {B<:StandardBasis,T,S,X} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})
function ⊗(p::ImmutableDensePolynomial{B,T,X,N},
           q::ImmutableDensePolynomial{B,S,X,M}) where {B<:StandardBasis,T,S,X,N,M}
    cs = fastconv(p.coeffs, q.coeffs)
    R = eltype(cs)
    ImmutableDensePolynomial{B,R,X,N+M-1}(cs)
end


#
function polynomial_composition(p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where {B<:StandardBasis,T,S,X,N,M}
    P = ImmutableDensePolynomial{B, promote_type(T,S), X, N*M}
    cs = evalpoly(q, p.coeffs)
    convert(P, cs)
end

# special cases of polynomial composition
# ... TBD ...

derivative(p::ImmutableDensePolynomial{B,T,X,0}) where {B<:StandardBasis,T,X} = p
function derivative(p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,N}
    N == 0 && return p
    hasnan(p) && return ⟒(p)(zero(T)/zero(T),X) # NaN{T}
    cs = ntuple(i -> i*p.coeffs[i+1], Val(N-1))
    R = eltype(cs)
    ImmutableDensePolynomial{StandardBasis,R,X,N-1}(cs)
end

integrate(p::ImmutableDensePolynomial{StandardBasis,T,X,0}) where {T,X} =
    ImmutableDensePolynomial{StandardBasis,Base.promote_op(/,T,Int),X,1}((0/1,))
function integrate(p::ImmutableDensePolynomial{StandardBasis,T,X,N}) where {T,X,N}
    N == 0 && return p # different type
    hasnan(p) && return ⟒(p)(zero(T)/zero(T), X) # NaN{T}
    z = zero(first(p.coeffs))
    cs = ntuple(i -> i > 1 ? p.coeffs[i-1]/(i-1) : z/1, Val(N+1))
    R = eltype(cs)
    ImmutableDensePolynomial{StandardBasis,R,X,N+1}(cs)
end
