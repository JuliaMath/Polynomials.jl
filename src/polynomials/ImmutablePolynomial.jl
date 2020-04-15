export ImmutablePolynomial

"""
    ImmutablePolynomial{T<:Number}(coeffs::AbstractVector{T}, var=:x)

Construct an immutable (static) polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct this through
`ImmutablePolynomial((a_0, a_1, ..., a_n))`.

The usual arithmetic operators are overloaded to work with polynomials as well as
with combinations of polynomials and scalars. However, operations involving two
polynomials of different variables causes an error.

This has the advantage over `Polynomial` as it can take advantage of faster polynomial evaluation
provided by `evalpoly` from Julia 1.4.

# Examples

```jldoctest
julia> ImmutablePolynomial([1, 0, 3, 4])
ImmutablePolynomial(1 + 3*x^2 + 4*x^3)

julia> ImmutablePolynomial([1, 2, 3], :s)
ImmutablePolynomial(1 + 2*s + 3*s^2)

julia> one(ImmutablePolynomial)
ImmutablePolynomial(1.0)
```

!!! note
    This was modeled after https://github.com/tkoolen/StaticUnivariatePolynomials.jl by @tkoolen.

"""
struct ImmutablePolynomial{N, T <: Number} <: StandardBasisPolynomial{T}
    coeffs::NTuple{N, T}
    var::Symbol
    function ImmutablePolynomial{N,T}(coeffs::NTuple{N,T}, var::Symbol=:x) where {N, T <: Number}
        new{N,T}(coeffs, var)
    end
    function ImmutablePolynomial{N, T}(coeffs::AbstractVector{S}, var::Symbol=:x) where {N, T <: Number, S}
        R = promote_type(T, S)
        length(coeffs) == 0 && return new{1, R}(NTuple(1, zero(R)), var)
        last_nz = findlast(!iszero, coeffs)
        M = max(1, last_nz === nothing ? 0 : last_nz)
        NN = max(N,M)
        if M < N
            cs = NTuple{NN,R}(i <= M ? R(coeffs[i]) : zero(R) for i in 1:NN)
        else
            cs = NTuple{NN,R}(R(c) for c in coeffs)
        end
        new{NN,R}(cs, var)
    end
end

## Can't do this so quickly, as we keep the parameter N
@register1 ImmutablePolynomial

# Base.convert(::Type{P}, p::P) where {P <: ImmutablePolynomial} = p
# function Base.convert(::Type{ImmutablePolynomial{N,T}}, p::ImmutablePolynomial{M,S}) where {N,T,M,S}
#     N >= M || throw(ArgumentError("XXX"))
#     ImmutablePolynomial{N,T}(NTuple{N,T}(T(p[i]) for i in 0:N-1), p.var)
# end
# Base.convert(P::Type{<:ImmutablePolynomial}, p::ImmutablePolynomial)  = p
# Base.convert(P::Type{<:ImmutablePolynomial}, p::Polynomial)  = P(coeffs(p), p.var)
# Base.convert(P::Type{<:Polynomial}, p::ImmutablePolynomial{N,T}) where {N,T}  = P(T[coeffs(p)...], p.var)
# Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{ImmutablePolynomial{N,S}}) where  {N, T,S} =
#     ImmutablePolynomial{N, promote_type(T, S)}
# Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{ImmutablePolynomial{M,S}}) where  {N, T,M,S} =
#     ImmutablePolynomial{max(N,M), promote_type(T, S)}
# Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{P}) where {N,T, P<:AbstractPolynomial} =
#     P
# Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{S}) where  {N,T,S<:Number}  =
#     ImmutablePolynomial{N, promote_type(T, S)}

# # various constructions based on  N, N&T, coefficients
# function ImmutablePolynomial{N,T}(coeffs::NTuple{M,S}, var::Polynomials.SymbolLike = :x) where{N,T, M,S} 
#     N >= M  || throw(ArgumentError(""))
#     ImmutablePolynomial{N,T}(NTuple{N, T}(i <= M ? T(coeffs[i]) : zero(T) for i in 1:N), Symbol(var))
# end
# ImmutablePolynomial{N}(coeffs::NTuple{M,T}, var::Polynomials.SymbolLike = :x) where{N,T, M} =
#     ImmutablePolynomial{M,T}(coeffs, Symbol(var))
# function ImmutablePolynomial{N}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where{N,T}
#     M = length(coeffs)
#     N >= M || throw(ArgumentError(""))
#     ImmutablePolynomial{N,T}(NTuple{N,T}(i <= M ? T(coeffs[i]) : zero(T) for i in  1:N), Symbol(var))
# end
ImmutablePolynomial(coeffs::NTuple{N,T}, var::Polynomials.SymbolLike = :x) where{N,T} =
     ImmutablePolynomial{N,T}(coeffs, Symbol(var))
function ImmutablePolynomial(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where{N,T}
    M = length(coeffs)
    ImmutablePolynomial{M,T}(NTuple{M,T}(x for x in coeffs), Symbol(var))
end

# # number
# ImmutablePolynomial{N,T}(n::Number, var::Polynomials.SymbolLike = :x) where {N,T} =
#     ImmutablePolynomial{N,T}(NTuple{N,T}(i==1 ? T(n) : zero(T) for i in 1:N), Symbol(var))
# ImmutablePolynomial{N}(n::T, var::Polynomials.SymbolLike = :x) where {N, T <: Number} =
#     ImmutablePolynomial{N,T}(NTuple{N,T}(i==1 ? n : zero(T) for i in 1:N), Symbol(var))
# ImmutablePolynomial(n::T, var::Polynomials.SymbolLike = :x) where {T <: Number} =
#     ImmutablePolynomial{1,T}(NTuple{1,T}(n))

# # P(var)  constructor
# ImmutablePolynomial{N,T}(var::Polynomials.SymbolLike = :x) where {N, T} = variable(ImmutablePolynomial{N,T}, var)
# ImmutablePolynomial{N}(var::Polynomials.SymbolLike = :x) where {N} = variable(ImmutablePolynomial{N}, var)
# ImmutablePolynomial(var::Polynomials.SymbolLike = :x)  = variable(ImmutablePolynomial{2, Int}, var)

# overrides from common.jl due to coeffs possibly being padded
Base.copy(p::P) where {P <: ImmutablePolynomial} = P(coeffs(p), p.var)
function Base.hash(p::ImmutablePolynomial{N,T}, h::UInt) where {N,T}
    n = findlast(!iszero, coeffs(p))
    isnothing(n) && return hash(p.var, hash(NTuple{0,T}(),h))
    hash(p.var, hash(coeffs(p)[1:n], h))
end

Base.one(::Type{ImmutablePolynomial{N,T}}) where {N,T}  = ImmutablePolynomial(tuple((i==1 ? one(T) : zero(T) for i in  1:N)...))
Base.zero(::Type{ImmutablePolynomial{N,T}}) where {N,T} = ImmutablePolynomial(zeros(T, N))

function degree(p::ImmutablePolynomial{N,T}) where {N, T}
    n = findlast(!iszero, coeffs(p))
    isnothing(n) ? -1 : n-1
end

function Base.:(==)(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    (p1.var == p2.var) || return false
    p1s, p2s = coeffs(p1), coeffs(p2)
    (N == M  && p1s == p2s) &&  return  true
    n1 = findlast(!iszero, p1s)
    n2 = findlast(!iszero, p2s)
    isnothing(n1) && isnothing(n2) && return true
    (isnothing(n1) || isnothing(n2)) && return false
    p1s[1:n1] == p2s[1:n2] &&  return true
    false
end
    
function Base.isapprox(p1::ImmutablePolynomial{N,T},
                       p2::ImmutablePolynomial{M,S};
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {N,T,M,S}
    p1.var == p2.var || error("p1 and p2 must have same var")
    NN = max(N,M)
    for i in 1:NN-1
        isapprox(p1[i],p2[i], rtol=rtol, atol=atol) || return false
    end
    true
end
function Base.isapprox(p1::ImmutablePolynomial{N,T}, n::S;
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {N, T,S}
    N >= 1 || return false
    isapprox(p1[0], n, rtol=rtol, atol=atol) || return false
    for i in 1:N-1
        isapprox(p1[i],zero(T), rtol=rtol, atol=atol) || return false
    end
    return true
end
function Base.chop(p::ImmutablePolynomial{N,T};
              rtol::Real = Base.rtoldefault(real(T)),
              atol::Real = 0)  where {N, T}
    cs = coeffs(p)
    for i in N:-1:1
        if !isapprox(cs[i], zero(T), rtol=rtol, atol=atol)
            return ImmutablePolynomial{i,T}(cs[1:i], p.var)
        end
    end
    zero(ImmutablePolynomial{0,T})
end
chop!(p::ImmutablePolynomial; kwargs...) =  chop(p; kwargs...)
function Base.truncate(p::ImmutablePolynomial{N,T};
                       rtol::Real = Base.rtoldefault(real(T)),
                       atol::Real = 0)  where {N, T}
    q = chop(p, rtol=rtol, atol=atol)
    iszero(q) && return  q
    cs = coeffs(q)
    thresh = maximum(abs,cs) * rtol + atol
    ImmutablePolynomial{length(cs), T}(map(c->abs(c) <= thresh ? zero(T) : c, coeffs(q)), p.var)
end
truncate!(p::ImmutablePolynomial; kwargs...) =  truncate(p; kwargs...)

LinearAlgebra.conj(p::P) where {P <: ImmutablePolynomial} = P(conj([aᵢ for aᵢ in coeffs(p)]))


##
## --------------------

(p::ImmutablePolynomial{N, T})(x::S) where {N, T,S} = evalpoly(x, coeffs(p))

function Base.:+(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    R = promote_type(T,S)
    m = max(N, M)
    return ImmutablePolynomial{m, R}(NTuple{m, R}(p1[i] + p2[i] for i in 0:m-1), p1.var)
end


function Base.:+(p::ImmutablePolynomial{N, T}, c::S) where {N, T,S<:Number}
    R = promote_type(T, S)
    cs = R[pi for pi in coeffs(p)]
    if isempty(cs)
        return ImmutablePolynomial{1,R}(NTuple{1,R}(R(c)), p.var)
    else
        cs[1] += c
    end
    ImmutablePolynomial{N,R}(NTuple{N,R}(c for c in cs), p.var)
end

## from https://github.com/tkoolen/StaticUnivariatePolynomials.jl/blob/master/src/monomial_basis.jl
@generated function Base.:*(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    P = M + N - 1
    exprs = Any[nothing for i = 1 : P]
    for i in 0 : M-1
        for j in 0 : N-1
            k = i + j + 1
            if exprs[k] === nothing
                exprs[k] = :(p1[$i] * p2[$j])
            else
                exprs[k] = :(muladd(p1[$i], p2[$j], $(exprs[k])))
            end
        end
    end
    return quote
        Base.@_inline_meta
        ImmutablePolynomial(tuple($(exprs...)))
    end
end



function Base.:*(p::ImmutablePolynomial{N,T}, c::S) where {N,T,S <: Number}
    R = Base.promote_op(*, T, S)
    return ImmutablePolynomial{N,R}(NTuple{N,R}(p[i]*c for i in eachindex(p)), p.var)
end

Base.:-(p::ImmutablePolynomial{N,T}) where {N,T} = ImmutablePolynomial(NTuple{N,T}(-pi for pi in p.coeffs), p.var)
Base.to_power_type(p::ImmutablePolynomial{N,T}) where {N,T} = p
