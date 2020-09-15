export ImmutablePolynomial

"""
    ImmutablePolynomial{T<:Number}(coeffs::AbstractVector{T}, var=:x)

Construct an immutable (static) polynomial from its coefficients `a`,
lowest order first, optionally in terms of the given variable `x`
where `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct
this through `ImmutablePolynomial((a_0, a_1, ..., a_n))` (assuming
`a_n ≠ 0`). As well, a vector or number can be used for construction.

The usual arithmetic operators are overloaded to work with polynomials
as well as with combinations of polynomials and scalars. However,
operations involving two polynomials of different variables causes an
error, though for `+` and `*` operations, constant polynomials are
treated as having no variable. 

As the coefficient size is a compile-time constant, immutable
polynomials can take advantage of faster polynomial evaluation
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
struct ImmutablePolynomial{T <: Number,  N} <: StandardBasisPolynomial{T}
    coeffs::NTuple{N, T}
    var::Symbol
    function ImmutablePolynomial{T,N}(coeffs::NTuple{N,T}, var::SymbolLike=:x) where {T <: Number,N}
        N == 0 && return new{T,0}(coeffs, Symbol(var))
        iszero(coeffs[end]) &&  throw(ArgumentError("Leading  term must  be  non-zero"))
        new{T,N}(coeffs, Symbol(var))
    end
    function ImmutablePolynomial{T,N}(coeffs::AbstractVector{S}, var::SymbolLike=:x) where {T <: Number, N, S}
        M = findlast(!iszero, coeffs)
        M == N || throw(ArgumentError("Leading  term must  be  non-zero; length of coefficients must be N"))
        cs = NTuple{N,T}(T(c) for c in coeffs[firstindex(coeffs):N])
        new{T,N}(cs, var)
    end
    function ImmutablePolynomial{T}(coeffs::NTuple{M,S}, var::SymbolLike=:x) where {T, S<: Number, M}
        N = findlast(!iszero, coeffs)
        if N == nothing
            return zero(ImmutablePolynomial{T}, var)
        else
            cs = NTuple{N,T}(coeffs[i] for  i in  1:N)
        end
        new{T,N}(cs, var)
    end
    function ImmutablePolynomial{T}(coeffs::Vector{M}, var::SymbolLike=:x) where {T, M}
        N = findlast(!iszero, coeffs)
        if N == nothing
            return zero(ImmutablePolynomial{T}, var)
        else
            cs = NTuple{N,T}(coeffs[i] for  i in  1:N)
        end
        new{T,N}(cs, var)
    end
    
end

@register ImmutablePolynomial

function ImmutablePolynomial{T,N}(coeffs::Tuple, var)  where {T,N}
    ImmutablePolynomial{T,N}(NTuple{N,T}(c for c in coeffs), var)
end


function ImmutablePolynomial(coeffs::NTuple{M,T}, var::SymbolLike=:x) where {M, T}
    N = findlast(!iszero, coeffs)
    if N == nothing
        return zero(ImmutablePolynomial, var) # no  type  passed in NTuple{M,T}
    else
        cs = NTuple{N,T}(T(c) for  c in  coeffs[1:N])
        end
    ImmutablePolynomial{T,N}(cs, var)
end

# Convenience; pass tuple to Polynomial
# Not documented, not sure this is a good idea as P(...)::P is not true...
Polynomial(coeffs::NTuple{N,T}, var::SymbolLike = :x) where{N,T} =
    ImmutablePolynomial(coeffs, var)
function Polynomial{T}(coeffs::NTuple{N,S}, var::SymbolLike = :x) where{N,T,S}
    ImmutablePolynomial{N,T}(T.(coeffs), var)
end

##
## ----
##
# overrides from common.jl due to coeffs possibly being padded, coeffs being non mutable, ...

Base.collect(p::P) where {P <: ImmutablePolynomial} = [pᵢ for pᵢ ∈ p]

Base.copy(p::P) where {P <: ImmutablePolynomial} = P(coeffs(p), p.var)

# catch q == 0 case
LinearAlgebra.norm(q::ImmutablePolynomial{T}, p::Real = 2) where {T} = degree(q) == -1 ? zero(T) : norm(coeffs(q), p)

#  zero, one, variable
function Base.zero(P::Type{<:ImmutablePolynomial},var::SymbolLike=:x)
    R = eltype(P)
    ImmutablePolynomial{R,0}(NTuple{0,R}(),var)
end

function  Base.one(P::Type{<:ImmutablePolynomial},var::SymbolLike=:x)
    R = eltype(P)
    ImmutablePolynomial{R,1}(NTuple{1,R}(1),var)
end
function variable(P::Type{<:ImmutablePolynomial},var::SymbolLike=:x)
    R  = eltype(P)
    ImmutablePolynomial{R,2}(NTuple{2,R}((0,1)),var)
end

# degree, isconstant
degree(p::ImmutablePolynomial{T,N}) where {T,N} = N - 1
isconstant(p::ImmutablePolynomial{T,N}) where {T,N}  = N <= 1

function Base.getindex(p::ImmutablePolynomial{T,N}, idx::Int) where {T <: Number,N}
    (idx <  0 || idx > N-1) && return zero(T)
    return p.coeffs[idx + 1]
end

Base.setindex!(p::ImmutablePolynomial, val::Number,  idx::Int) = throw(ArgumentError("Immutable polynomials are immutable"))

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

##
## Change ≈ to handle tuples for coefficients *and* get handling of Inf correct
## right now
## Inf ≈ Inf # true
## [Inf] ≈ [Inf] # true
## P([Inf]) ≈ P([Inf]) # false
## P([Inf]) ≈ Inf # false
## This fixes the last two cases for P=ImmutablePolynomial, and could replace
## isapprox in common.jl
##
## check Inf values (matching Vector ≈ Vector)
## @test P([Inf]) ≈ Inf
## @test !(P([Inf]) ≈ P([-Inf])) # default compares zero(P) to zero(P)
## @test !(P([1,2,Inf]) ≈ P([1,3,Inf])) # default uses truncate with  Inf which clobbers numbers
function Base.isapprox(p1::ImmutablePolynomial{T},
    p2::ImmutablePolynomial{S};
    rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {T,S}
    
    p1, p2 = promote(p1, p2)
    check_same_variable(p1, p2)  || error("p1 and p2 must have same var")

    # copy over from abstractarray.jl
    Δ  = norm(p1-p2)
    if isfinite(Δ)
        return Δ <= max(atol, rtol*max(norm(p1), norm(p2)))
    else
        for i in 0:max(degree(p1), degree(p2))
            isapprox(p1[i], p2[i]; rtol=rtol, atol=atol) || return false
        end
        return true
    end
end

function Base.isapprox(p1::P,
                       n::S;
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {T,S, P<:ImmutablePolynomial{T}}
    return isapprox(p1, ⟒(P){T}(n,p1.var))
end

function Base.chop(p::ImmutablePolynomial{T,N};
              rtol::Real = Base.rtoldefault(real(T)),
              atol::Real = 0)  where {T,N}
    cs = coeffs(p)
    for i in N:-1:1
        if !isapprox(cs[i], zero(T), rtol=rtol, atol=atol)
            return ImmutablePolynomial{T,i}(cs[1:i], p.var)
        end
    end
    zero(ImmutablePolynomial{T})
end

function Base.truncate(p::ImmutablePolynomial{T,N};
                       rtol::Real = Base.rtoldefault(real(T)),
                       atol::Real = 0)  where {T,N}
    q = chop(p, rtol=rtol, atol=atol)
    iszero(q) && return  q
    cs = coeffs(q)
    thresh = maximum(abs,cs) * rtol + atol
    ImmutablePolynomial(map(c->abs(c) <= thresh ? zero(T) : c, coeffs(q)), p.var)
end

# no in-place chop! and truncate!
chop!(p::ImmutablePolynomial; kwargs...) =  chop(p; kwargs...)
truncate!(p::ImmutablePolynomial; kwargs...) =  truncate(p; kwargs...)

LinearAlgebra.conj(p::P) where {P <: ImmutablePolynomial} = P(conj([aᵢ for aᵢ in coeffs(p)]))


##
## --------------------
##

(p::ImmutablePolynomial{T,N})(x::S) where {T,N,S} = evalpoly(x, coeffs(p))


function Base.:+(p1::ImmutablePolynomial{T,N}, p2::ImmutablePolynomial{S,M}) where {T,N,S,M}

    R = promote_type(S,T)
    iszero(N) && return ImmutablePolynomial{R}(coeffs(p2), p2.var)
    iszero(M) && return ImmutablePolynomial{R}(coeffs(p1), p1.var)
    
    isconstant(p1) && p1.var != p2.var && return p2 + p1[0]*one(ImmutablePolynomial{T},p2.var)
    isconstant(p2) && p1.var != p2.var && return p1 + p2[0]*one(ImmutablePolynomial{S},p1.var)

    p1.var != p2.var && error("Polynomials must have same variable")


    if  N == M
        cs = coeffs(p1) .+ coeffs(p2)
        return ImmutablePolynomial(cs, p1.var) # N  unknown, as leading terms can cancel
    elseif N < M
        ⊕(p2,p1)
    else
        ⊕(p1,p2)
    end

end

# assume N > M.
# As N ≠ M, we are assured of size of output (max(N,M)), so we generate the function
@generated function ⊕(p1::ImmutablePolynomial{T,N}, p2::ImmutablePolynomial{S,M}) where {T,N,S,M}
    R = promote_type(T,S)

    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(p1.coeffs[$i] + p2.coeffs[$i])
    end
    for i in M+1:N
        exprs[i] =:(p1.coeffs[$i])
    end

    return quote
        Base.@_inline_meta
        ImmutablePolynomial{$(R),$(N)}(tuple($(exprs...)), p1.var)
    end

end



function Base.:*(p1::ImmutablePolynomial, p2::ImmutablePolynomial)
    isconstant(p1) && return p2 * p1[0] 
    isconstant(p2) && return p1 * p2[0]
    p1.var != p2.var && error("Polynomials must have same variable")
    p1 ⊗ p2
end

## Static size of product makes generated functions  a good choice
## from https://github.com/tkoolen/StaticUnivariatePolynomials.jl/blob/master/src/monomial_basis.jl
@generated function ⊗(p1::ImmutablePolynomial{T,N}, p2::ImmutablePolynomial{S,M}) where {T,N,S,M}
    P = M + N - 1
    R = promote_type(T,S)
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
        ImmutablePolynomial{$(R),$(max(0,P))}(tuple($(exprs...)), p1.var)
    end

end

# scalar ops
function Base.:+(p::ImmutablePolynomial{T,N}, c::S) where {T, N, S<:Number}
    R = promote_type(T,S)
    iszero(c) && return ImmutablePolynomial{R,N}(R.(p.coeffs), p.var)
    N == 0 && return ImmutablePolynomial{R,1}((R(c),), p.var)
    N == 1 && return ImmutablePolynomial(R[p[0]+c])

    return p + ImmutablePolynomial{R,1}(R[c],p.var)
    
    cs = NTuple{N,R}(iszero(i) ? p[i]+c : p[i] for i in 0:N-1)
    return ImmutablePolynomial{R,N}(cs, p.var)
end

function Base.:*(p::ImmutablePolynomial{T,N}, c::S) where {T, N, S <: Number}
    R = promote_type(T,S)
    iszero(c)  && return zero(ImmutablePolynomial{R})
    ImmutablePolynomial{R,N}(p.coeffs .* c, p.var)
end

function Base.:/(p::ImmutablePolynomial{T,N}, c::S) where {T,N,S <: Number}
    R = eltype(one(T)/one(S))
    isinf(c)  && return zero(ImmutablePolynomial{R})
    ImmutablePolynomial{R,N}(p.coeffs ./ c, p.var)
end

Base.:-(p::ImmutablePolynomial{T,N}) where {T,N} = ImmutablePolynomial(NTuple{N,T}(-pi for pi in p.coeffs), p.var)

Base.to_power_type(p::ImmutablePolynomial{T,N}) where {T,N} = p
