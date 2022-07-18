export ImmutablePolynomial

"""
    ImmutablePolynomial{T, X, N}(coeffs::AbstractVector{T})

Construct an immutable (static) polynomial from its coefficients
`a₀, a₁, …, aₙ`,
lowest order first, optionally in terms of the given variable `x`
where `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct
this through `ImmutablePolynomial((a_0, a_1, ..., a_n))` (assuming
`a_n ≠ 0`). As well, a vector or number can be used for construction.


The usual arithmetic operators are overloaded to work with polynomials
as well as with combinations of polynomials and scalars. However,
operations involving two non-constant polynomials of different variables causes an
error. Unlike other polynomials, `setindex!` is not defined for `ImmutablePolynomials`.

As the degree of the polynomial (`+1`) is a compile-time constant,
several performance improvements are possible. For example, immutable
polynomials can take advantage of faster polynomial evaluation
provided by `evalpoly` from Julia 1.4; similar methods are also used
for addtion and multiplication.

However, as the degree is included in the type, promotion between
immutable polynomials can not promote to a common type. As such, they
are precluded from use in rational functions.

!!! note
    `ImmutablePolynomial` is not axis-aware, and it treats `coeffs` simply as a list of coefficients with the first
    index always corresponding to the constant term.

# Examples

```jldoctest
julia> using  Polynomials

julia> ImmutablePolynomial((1, 0, 3, 4))
ImmutablePolynomial(1 + 3*x^2 + 4*x^3)

julia> ImmutablePolynomial((1, 2, 3), :s)
ImmutablePolynomial(1 + 2*s + 3*s^2)

julia> one(ImmutablePolynomial)
ImmutablePolynomial(1.0)
```

!!! note
    This was modeled after [StaticUnivariatePolynomials](https://github.com/tkoolen/StaticUnivariatePolynomials.jl) by `@tkoolen`.

"""
struct ImmutablePolynomial{T, X, N} <: StandardBasisPolynomial{T, X}
    coeffs::NTuple{N, T}
    function ImmutablePolynomial{T,X,N}(coeffs::NTuple{N,T}) where {T, X, N}
        N == 0 && return new{T,X, 0}(coeffs)
        iszero(coeffs[end]) &&  throw(ArgumentError("Leading  term must  be  non-zero"))
        new{T,X,N}(coeffs)
    end
end

@register ImmutablePolynomial

## Various interfaces
## Abstract Vector coefficients
function ImmutablePolynomial{T,X, N}(coeffs::AbstractVector{T})  where {T,X, N}
    cs = NTuple{N,T}(coeffs[i] for i ∈ firstindex(coeffs):N)
    ImmutablePolynomial{T, X, N}(cs)

end

function ImmutablePolynomial{T,X, N}(coeffs::AbstractVector{S})  where {T,X, N, S}
    cs = NTuple{N,T}(coeffs[i] for i ∈ firstindex(coeffs):N)
    ImmutablePolynomial{T, X, N}(cs)

end

function ImmutablePolynomial{T,X}(coeffs::AbstractVector{S})  where {T,X,S}
    R = promote_type(T,S)

    if Base.has_offset_axes(coeffs)
        @warn "ignoring the axis offset of the coefficient vector"
    end
    N = findlast(!iszero, coeffs)
    N == nothing && return ImmutablePolynomial{R,X,0}(())
    N′ = N + 1 - firstindex(coeffs)
    ImmutablePolynomial{T, X, N′}([coeffs[i] for i ∈ firstindex(coeffs):N])
end

## -- Tuple arguments
function ImmutablePolynomial{T,X}(coeffs::Tuple)  where {T,X}
    N = findlast(!iszero, coeffs)
    N == nothing && return zero(ImmutablePolynomial{T,X})
    ImmutablePolynomial{T,X,N}(NTuple{N,T}(coeffs[i] for i in 1:N))
end

ImmutablePolynomial{T}(coeffs::Tuple, var::SymbolLike=:x)  where {T} = ImmutablePolynomial{T,Symbol(var)}(coeffs)

function ImmutablePolynomial(coeffs::Tuple, var::SymbolLike=:x)
    cs = NTuple(promote(coeffs...))
    T = eltype(cs)
    ImmutablePolynomial{T, Symbol(var)}(cs)
end

##
## ----
##
# overrides from common.jl due to  coeffs being non mutable, N in type parameters

Base.copy(p::P) where {P <: ImmutablePolynomial} = P(coeffs(p))
Base.similar(p::ImmutablePolynomial, args...) =
    similar(collect(oeffs(p)), args...)
# degree, isconstant
degree(p::ImmutablePolynomial{T,X, N}) where {T,X,N} = N - 1 # no trailing zeros
isconstant(p::ImmutablePolynomial{T,X,N}) where {T,X,N}  = N <= 1

Base.setindex!(p::ImmutablePolynomial, val,  idx::Int) = throw(ArgumentError("ImmutablePolynomials are immutable"))

for op in [:isequal, :(==)]
    @eval function Base.$op(p1::ImmutablePolynomial{T,N}, p2::ImmutablePolynomial{S,M}) where {T,N,S,M}
        check_same_variable(p1,p2) || return false
        p1s, p2s = coeffs(p1), coeffs(p2)
        (N == M  && $op(p1s,p2s)) &&  return  true
        n1 = findlast(!iszero, p1s) # now trim out zeros
        n2 = findlast(!iszero, p2s)
        (n1 == nothing && n2 == nothing) && return true
        (n1 == nothing || n2  == nothing) && return false
        $op(p1s[1:n1],p2s[1:n2]) &&  return true
        false
    end
end

# in common.jl these call chop! and truncate!
function Base.chop(p::ImmutablePolynomial{T,X};
              rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0)  where {T,X}
    ps = chop(p.coeffs; rtol=rtol, atol=atol)
    return ImmutablePolynomial{T,X}(ps)
end

function Base.truncate(p::ImmutablePolynomial{T,X};
                  rtol::Real = Base.rtoldefault(real(T)),
                  atol::Real = 0)  where {T,X}
    ps = truncate(p.coeffs; rtol=rtol, atol=atol)
    ImmutablePolynomial{T,X}(ps)
end


##
## --------------------
##

## Addition
# scalar ops
function Base.:+(p::P, c::S) where {T, X, N, P <: ImmutablePolynomial{T,X,N}, S<:Number}
    R = promote_type(T,S)

    iszero(c) && return ImmutablePolynomial{R,X,N}(convert(NTuple{N,R},p.coeffs))
    N == 0 && return ImmutablePolynomial{R,X,1}(NTuple{1,R}(c))
    N == 1 && return ImmutablePolynomial((p[0]+c,), X)

    cs = ⊕(P, convert(NTuple{N,R},p.coeffs), NTuple{1,R}(c))
    q = ImmutablePolynomial{R,X,N}(cs)

    return q

end

Base.:-(p::ImmutablePolynomial{T,X,N}) where {T,X,N} = ImmutablePolynomial{T,X,N}(.-p.coeffs)

function Base.:+(p1::P, p2::Q) where {T,X,N,P<:ImmutablePolynomial{T,X,N},
                                      S,  M,Q<:ImmutablePolynomial{S,X,M}}

    R = promote_type(T,S)
    P′ = ImmutablePolynomial{R,X}
    if  N == M
        cs = ⊕(P, p1.coeffs, p2.coeffs)
        return P′(R.(cs))
    elseif N < M
        cs = ⊕(P, p2.coeffs, p1.coeffs)
        return P′{M}(R.(cs))
    else
        cs = ⊕(P, p1.coeffs, p2.coeffs)
        return P′{N}(R.(cs))
    end

end

## multiplication

function scalar_mult(p::ImmutablePolynomial{T,X,N}, c::S) where {T, X,N, S <: Number}
    R = eltype(p[0] * c * 0)
    (N == 0  || iszero(c)) && return zero(ImmutablePolynomial{R,X})
    cs = p.coeffs .* c
    return ImmutablePolynomial(cs, X)
 end

function scalar_mult(c::S, p::ImmutablePolynomial{T,X,N}) where {T, X,N, S <: Number}
    R = eltype(p[0] * c * 0)
    (N == 0  || iszero(c)) && return zero(ImmutablePolynomial{R,X})
    cs = p.coeffs .* c
    return ImmutablePolynomial(cs, X)
end

function Base.:/(p::ImmutablePolynomial{T,X,N}, c::S) where {T,X,N,S<:Number}
    R = eltype(one(T)/one(S))
    P = ImmutablePolynomial{R,X}
    (N == 0  || isinf(c)) && return zero(P)
    cs = p.coeffs ./ c
    iszero(cs[end]) ? P(cs) : P{N}(cs) # more performant to specify when N is known
end


function Base.:*(p1::ImmutablePolynomial{T,X,N}, p2::ImmutablePolynomial{S,X,M}) where {T,S,X,N,M}
    R = promote_type(T,S)
    P = ImmutablePolynomial{R,X}

    (iszero(N) || iszero(M)) && return zero(P)

    cs = ⊗(ImmutablePolynomial, p1.coeffs, p2.coeffs) #(p1.coeffs) ⊗ (p2.coeffs)
    iszero(cs[end]) ? P(cs) : P{N+M-1}(cs)  # more performant to specify when N is known

end

Base.to_power_type(p::ImmutablePolynomial{T,X,N}) where {T,X,N} = p


## more performant versions borrowed from StaticArrays
## https://github.com/JuliaArrays/StaticArrays.jl/blob/master/src/linalg.jl
LinearAlgebra.norm(q::ImmutablePolynomial{T,X,0}) where {T,X} = zero(real(float(T)))
LinearAlgebra.norm(q::ImmutablePolynomial) = _norm(q.coeffs)
LinearAlgebra.norm(q::ImmutablePolynomial, p::Real) = _norm(q.coeffs, p)

@generated function _norm(a::NTuple{N,T}) where {T, N}

    expr = :(abs2(a[1]))
    for j = 2:N
        expr = :($expr + abs2(a[$j]))
    end

    return quote
        $(Expr(:meta, :inline)) # 1.8 deprecation
        #Base.@inline
        @inbounds return sqrt($expr)
    end

end

_norm_p0(x) = iszero(x) ? zero(x) : one(x)
@generated function _norm(a::NTuple{N,T}, p::Real) where {T, N}
    expr = :(abs(a[1])^p)
    for j = 2:N
        expr = :($expr + abs(a[$j])^p)
    end

    expr_p1 = :(abs(a[1]))
    for j = 2:N
        expr_p1 = :($expr_p1 + abs(a[$j]))
    end

    return quote
        $(Expr(:meta, :inline)) # 1.8 deprecation
        #Base.@inline
        if p == Inf
            return mapreduce(abs, max, a)
        elseif p == 1
            @inbounds return $expr_p1
        elseif p == 2
            return norm(a)
        elseif p == 0
            return mapreduce(_norm_p0, +, a)
        else
            @inbounds return ($expr)^(inv(p))
        end
    end

end
