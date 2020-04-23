export ImmutablePolynomial

"""
    ImmutablePolynomial{T<:Number}(coeffs::AbstractVector{T}, var=:x)

Construct an immutable (static) polynomial from its coefficients `a`,
lowest order first, optionally in terms of the given variable `x`
where `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct this through
`ImmutablePolynomial((a_0, a_1, ..., a_n))`.

The usual arithmetic operators are overloaded to work with polynomials
as well as with combinations of polynomials and scalars. However,
operations involving two polynomials of different variables causes an
error, though for `+` and `*` operations, constant polynomials are
treated as having no variable. (This adds a runtime check, but is useful for using with
matrices.)

This has the advantage over `Polynomial` as it can take advantage of faster polynomial evaluation
provided by `evalpoly` from Julia 1.4. 

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
    This was modeled after https://github.com/tkoolen/StaticUnivariatePolynomials.jl by @tkoolen.

"""
struct ImmutablePolynomial{N, T <: Number} <: StandardBasisPolynomial{T}
    coeffs::NTuple{N, T}
    var::Symbol
    function ImmutablePolynomial{N,T}(coeffs::NTuple{N,T}, var::Symbol=:x) where {N, T <: Number}
        new{N,T}(coeffs, var)
    end
    function ImmutablePolynomial{N,T}(coeffs::NTuple{M,S}, var::Symbol=:x) where {N, M, T, S<: Number}
        M > N && throw(ArgumentError("coeffs too big for the specified type"))
        if M < N
            cs = NTuple{N,T}(i <= M ? T(coeffs[i]) : zero(T) for i in 1:N)
        else
            cs = NTuple{N,T}(T(c) for c in coeffs)
        end
        new{N,T}(cs, var)
    end
    function ImmutablePolynomial{N,T}(coeffs::AbstractVector{S}, var::Symbol=:x) where {N, T <: Number, S}
        M = findlast(!iszero, coeffs)
        M == nothing && return zero(ImmutablePolynomial{N,T})
        if M < N
            cs = NTuple{N,T}(i <= M ? T(coeffs[i]) : zero(T) for i in 1:N)
        else
            cs = NTuple{N,T}(T(c) for c in coeffs)
        end
        new{N,T}(cs, var)
    end
end

@register1 ImmutablePolynomial

## promote N,M case
function Base.promote(p::ImmutablePolynomial{N,T}, q::ImmutablePolynomial{M,S}) where {N,T,M,S}
    NN, R = max(N,M), promote_type(T,S)
    ImmutablePolynomial{NN,R}(p.coeffs, p.var), ImmutablePolynomial{NN,R}(q.coeffs, q.var)
end
## Need defaults when N not specified
# tuple
ImmutablePolynomial(coeffs::NTuple{N,T}, var::Polynomials.SymbolLike = :x) where{N,T} =
     ImmutablePolynomial{N,T}(coeffs, Symbol(var))

# vector
function ImmutablePolynomial(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where{N,T}
    M = length(coeffs)
    M == 0 && return zero(ImmutablePolynomial{1, T}, var)
    ImmutablePolynomial{M,T}(NTuple{M,T}(x for x in coeffs), Symbol(var))
end

# number
ImmutablePolynomial(n::T, var::Polynomials.SymbolLike = :x) where {T <: Number} =
     ImmutablePolynomial{1,T}(NTuple{1,T}(n), Symbol(var))

# variable
ImmutablePolynomial(var::Polynomials.SymbolLike = :x)  = variable(ImmutablePolynomial{2, Int}, var)

# Convenience; pass tuple to Polynomial
# Not documented, not sure this is a good idea as P(...)::P is not true...
Polynomial(coeffs::NTuple{N,T}, var::Polynomials.SymbolLike = :x) where{N,T} =
    ImmutablePolynomial(coeffs, var)
function Polynomial{T}(coeffs::NTuple{N,S}, var::Polynomials.SymbolLike = :x) where{N,T,S}
    ImmutablePolynomial{N,T}(T.(coeffs), var)
end

##
## ----
##

# overrides from common.jl due to coeffs possibly being padded, coeffs being no mutable, ...
Base.promote_rule(::Type{<:ImmutablePolynomial{N, T}}, ::Type{<:ImmutablePolynomial{M,S}}) where {N,T,M,S} =
    ImmutablePolynomial{max(N,M), promote_type(T, S)}

Base.convert(::Type{<:ImmutablePolynomial{N, T}}, p::ImmutablePolynomial{M,S}) where {N,T,M,S} = 
    ImmutablePolynomial{N,T}(p.coeffs,  p.var)


Base.copy(p::P) where {P <: ImmutablePolynomial} = P(coeffs(p), p.var)

function Base.hash(p::ImmutablePolynomial{N,T}, h::UInt) where {N,T}
    n = findlast(!iszero, coeffs(p))
    n == nothing && return hash(p.var, hash(NTuple{0,T}(),h))
    hash(p.var, hash(coeffs(p)[1:n], h))
end

Base.one(::Type{ImmutablePolynomial{N,T}},
         var=:x) where {N,T}  = ImmutablePolynomial(tuple((i==1 ? one(T) : zero(T) for i in  1:N)...), var)
Base.one(::Type{ImmutablePolynomial{N}}, var=:x) where {N}  = one(ImmutablePolynomial{N,Int}, var)
Base.zero(::Type{ImmutablePolynomial{N,T}}, var=:x) where {N,T} = ImmutablePolynomial(zeros(T, N), var)
Base.zero(::Type{ImmutablePolynomial{N}}, var=:x) where {N}  = zero(ImmutablePolynomial{N,Int}, var)

function degree(p::ImmutablePolynomial{N,T}) where {N, T}
    n = findlast(!iszero, coeffs(p))
    n == nothing ? -1 : n-1
end

for op in [:isequal, :(==)]
    @eval function Base.$op(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
        (p1.var == p2.var) || return false
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
                       atol::Real = 0,) where {N, T,S <: Number}
    return isapprox(p1, ImmutablePolynomial(n, p1.var))
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
##

(p::ImmutablePolynomial{N, T})(x::S) where {N, T,S} = evalpoly(x, coeffs(p))

# used to treat constants as havinig same variable as counterpart in + and *
function promote_constants(p::ImmutablePolynomial{N,T}, q::ImmutablePolynomial{M,S}) where {N,T,M,S}
    if  degree(p) <= 0
        p  = ImmutablePolynomial{N,T}(p.coeffs, q.var)
    elseif degree(q) <= 0
        q  = ImmutablePolynomial{M,S}(q.coeffs, p.var)
    end
    p,q
end

function Base.:+(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    p1,p2 = promote_constants(p1, p2)
    p1.var != p2.var && error("Polynomials must have same variable")
    
    R = Base.promote_op(+, T,S)
    NN = max(N, M)
    P = ImmutablePolynomial{NN,R}
    ⊕(convert(P,p1), convert(P,p2))
end

function ⊕(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{N,T}) where {N,T}
    cs = NTuple{N, T}(p1[i] + p2[i] for i in 0:N-1)
    return ImmutablePolynomial{N, T}(cs, p1.var)
end


Base.:+(p::ImmutablePolynomial{N, T}, c::S) where {N, T,S<:Number} =
    p + ImmutablePolynomial((c,), p.var)

function Base.:*(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    p1,p2 = promote_constants(p1, p2)
    p1.var != p2.var && error("Polynomials must have same variable")
    p1 ⊗ p2
end

## from https://github.com/tkoolen/StaticUnivariatePolynomials.jl/blob/master/src/monomial_basis.jl
@generated function ⊗(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    R = eltype(one(T)*one(S))
    P = M + N - 1
    exprs = Any[nothing for i = 1 : P]
    for i in 1 : N
        for j in 1 : M
            k = i + j - 1
            if exprs[k] === nothing
                exprs[k] = :(p1.coeffs[$i] * p2.coeffs[$j])
            else
                exprs[k] = :(muladd(p1.coeffs[$i], p2.coeffs[$j], $(exprs[k])))
            end
        end
    end
    return quote
        Base.@_inline_meta
        ImmutablePolynomial{$(max(0,P)), $(R)}(tuple($(exprs...)), p1.var)
    end
end


function Base.:*(p::ImmutablePolynomial{N,T}, c::S) where {N,T,S <: Number}
    R = eltype(one(T)*one(S))
    return ImmutablePolynomial{N,R}(NTuple{N,R}(p[i]*c for i in eachindex(p)), p.var)
end

Base.:-(p::ImmutablePolynomial{N,T}) where {N,T} = ImmutablePolynomial(NTuple{N,T}(-pi for pi in p.coeffs), p.var)

Base.to_power_type(p::ImmutablePolynomial{N,T}) where {N,T} = p
