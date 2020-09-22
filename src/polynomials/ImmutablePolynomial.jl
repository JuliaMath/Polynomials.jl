export ImmutablePolynomial

"""
    ImmutablePolynomial{T<:Number, N}(coeffs::AbstractVector{T}, var=:x)

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

As the coefficient size is a compile-time constant, several performance
improvements are possible. For example, immutable polynomials can take advantage of 
faster polynomial evaluation provided by `evalpoly` from Julia 1.4.


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

    
end

@register ImmutablePolynomial

## Various interfaces
function ImmutablePolynomial{T,N}(coeffs::Tuple, var::SymbolLike=:x)  where {T,N}
    ImmutablePolynomial{T,N}(NTuple{N,T}(T.(coeffs)), var)
end

function ImmutablePolynomial{T,N}(coeffs::AbstractVector{S}, var::SymbolLike=:x) where {T <: Number, N, S}
    ImmutablePolynomial{T,N}(NTuple{N,T}(tuple(coeffs...)), var)
end

## --
function ImmutablePolynomial{T}(coeffs::NTuple{M,S}, var::SymbolLike=:x) where {T, S<: Number, M}
    N = findlast(!iszero, coeffs)
    if N == nothing
        return zero(ImmutablePolynomial{T}, var)
    else
        cs = NTuple{N,T}(coeffs[i] for  i in  1:N)
        return ImmutablePolynomial{T,N}(cs, var)
    end
end

function ImmutablePolynomial{T}(coeffs::Tuple, var::SymbolLike=:x)  where {T}
    ImmutablePolynomial{T}(T.(coeffs), var)
end

# entry point from abstract.jl; note T <: Number
function ImmutablePolynomial{T}(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {T <: Number}
    M = length(coeffs)
    ImmutablePolynomial{T}(NTuple{M,T}(tuple(coeffs...)), var)
end

function ImmutablePolynomial{T}(coeffs::OffsetArray{T,1, Array{T, 1}}, var::SymbolLike=:x) where {T <: Number}
    cs = zeros(T, 1 + lastindex(coeffs))
    cs[1 .+ (firstindex(coeffs):lastindex(coeffs))] = coeffs.parent
    ImmutablePolynomial{T}(cs, var)
end

## --

function ImmutablePolynomial(coeffs::Tuple, var::SymbolLike=:x)
    cs = NTuple(promote(coeffs...))
    T = eltype(cs)
    ImmutablePolynomial{T}(cs, var)
end

# # need to catch map case wihch return Variables, not Values
# function ImmutablePolynomial(coeffs::Cs, var::SymbolLike=:x) where {
#     M, T, Cs <: Union{Values{M,T}, Variables{M,T}}}
#     ImmutablePolynomial{T}(Values(coeffs.v), var)
# end




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
# overrides from common.jl due to  coeffs being non mutable, N in type parameters
Base.collect(p::P) where {P <: ImmutablePolynomial} = [pᵢ for pᵢ ∈ p]

Base.copy(p::P) where {P <: ImmutablePolynomial} = P(coeffs(p), p.var)

## defining these speeds things up
function Base.zero(P::Type{<:ImmutablePolynomial}, var::SymbolLike=:x)
    R = eltype(P)
    ImmutablePolynomial{R,0}(NTuple{0,R}(),var)
end

function  Base.one(P::Type{<:ImmutablePolynomial}, var::SymbolLike=:x)
    R = eltype(P)
    ImmutablePolynomial{R,1}(NTuple{1,R}(1),var)
end
function variable(P::Type{<:ImmutablePolynomial}, var::SymbolLike=:x)
    R  = eltype(P)
    ImmutablePolynomial{R,2}(NTuple{2,R}((0,1)),var)
end


# degree, isconstant
degree(p::ImmutablePolynomial{T,N}) where {T,N} = N - 1 # no trailing zeros
isconstant(p::ImmutablePolynomial{T,N}) where {T,N}  = N <= 1

function Base.getindex(p::ImmutablePolynomial{T,N}, idx::Int) where {T <: Number,N}
    (idx <  0 || idx > N-1) && return zero(T)
    return p.coeffs[idx + 1]
end

Base.setindex!(p::ImmutablePolynomial, val::Number,  idx::Int) = throw(ArgumentError("ImmutablePolynomials are immutable"))

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
function Base.chop(p::ImmutablePolynomial{T,N};
              rtol::Real = Base.rtoldefault(real(T)),
              atol::Real = 0)  where {T,N}
    cs = coeffs(p)
    for i in N:-1:1
        if !isapprox(cs[i], zero(T), rtol=rtol, atol=atol)
            return ImmutablePolynomial{T,i}(cs[1:i], p.var)
        end
    end
    zero(ImmutablePolynomial{T}, p.var)
end

function Base.truncate(p::ImmutablePolynomial{T,N};
                       rtol::Real = Base.rtoldefault(real(T)),
                       atol::Real = 0)  where {T,N}
    q = chop(p, rtol=rtol, atol=atol)
    iszero(q) && return  q
    cs = coeffs(q)
    thresh = maximum(abs,cs) * rtol + atol
    cs′ = map(c->abs(c) <= thresh ? zero(T) : c, cs)
    ImmutablePolynomial{T}(tuple(cs′...), p.var)
end

# no in-place chop! and truncate!
chop!(p::ImmutablePolynomial; kwargs...) =  chop(p; kwargs...)
truncate!(p::ImmutablePolynomial; kwargs...) =  truncate(p; kwargs...)

##
## --------------------
##

(p::ImmutablePolynomial{T,N})(x::S) where {T,N,S} = evalpoly(x, p.coeffs)


function Base.:+(p1::ImmutablePolynomial{T,N}, p2::ImmutablePolynomial{S,M}) where {T,N,S,M}

    R = promote_type(S,T)
    iszero(N) && return ImmutablePolynomial{R}(coeffs(p2), p2.var)
    iszero(M) && return ImmutablePolynomial{R}(coeffs(p1), p1.var)
    
    isconstant(p1) && p1.var != p2.var && return p2 + p1[0]*one(ImmutablePolynomial{R}, p2.var)
    isconstant(p2) && p1.var != p2.var && return p1 + p2[0]*one(ImmutablePolynomial{R}, p1.var)

    p1.var != p2.var && error("Polynomials must have same variable")


    if  N == M
        cs = NTuple{N,R}(p1[i] + p2[i] for i in 0:N-1)
        ImmutablePolynomial{R}(cs, p1.var)        
    elseif N < M
        cs = (p2.coeffs) ⊕ (p1.coeffs)
        ImmutablePolynomial{R,M}(cs, p1.var)        
    else
        cs = (p1.coeffs) ⊕ (p2.coeffs)
        ImmutablePolynomial{R,N}(cs, p1.var)                
    end

    

end


function Base.:*(p1::ImmutablePolynomial{T,N}, p2::ImmutablePolynomial{S,M}) where {T,N,S,M}
    isconstant(p1) && return p2 * p1[0] 
    isconstant(p2) && return p1 * p2[0]
    p1.var != p2.var && error("Polynomials must have same variable")
    R = promote_type(S,T)
    cs = (p1.coeffs) ⊗ (p2.coeffs)
    ImmutablePolynomial{R, N+M-1}(cs, p1.var)
end

# Padded vector sum of two tuples assuming N > M
# assume N > M.
# As N ≠ M, we are assured of size of output (max(N,M)), so we generate the function
@generated function ⊕(p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}

    R = promote_type(T,S)

    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(p1[$i] + p2[$i])
    end
    for i in M+1:N
        exprs[i] =:(p1[$i])
    end

    return quote
        Base.@_inline_meta
        tuple($(exprs...))
    end

end



## Static size of product makes generated functions  a good choice
## from https://github.com/tkoolen/StaticUnivariatePolynomials.jl/blob/master/src/monomial_basis.jl
## convolution of two tuples
@generated function ⊗(p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}
    P = M + N - 1
    R = promote_type(T,S)
    exprs = Any[nothing for i = 1 : P]
    for i in 1 : N
        for j in 1 : M
            k = i + j - 1
            if exprs[k] === nothing
                exprs[k] = :(p1[$i] * p2[$j])
            else
                exprs[k] = :(muladd(p1[$i], p2[$j], $(exprs[k])))
            end
        end
    end

    return quote
        Base.@_inline_meta
        tuple($(exprs...))        
    end

end

# scalar ops
function Base.:+(p::ImmutablePolynomial{T,N}, c::S) where {T, N, S<:Number}
    R = promote_type(T,S)

    iszero(c) && return ImmutablePolynomial{R,N}(p.coeffs, p.var)
    N == 0 && return ImmutablePolynomial{R,1}((c,), p.var)
    N == 1 && return ImmutablePolynomial((p[0]+c,), p.var)

#    cs = NTuple{N,R}(iszero(i) ? p[i]+c : p[i] for i in 0:N-1)
#    return ImmutablePolynomial{R,N}(cs, p.var)
    
    q = ImmutablePolynomial{R,1}((c,), p.var)
    return p + q

end

function Base.:*(p::ImmutablePolynomial{T,N}, c::S) where {T, N, S <: Number}
    R = promote_type(T,S)
    iszero(c) && return zero(ImmutablePolynomial{R}, p.var)
    ImmutablePolynomial{R,N}(p.coeffs .* c, p.var)
end

function Base.:/(p::ImmutablePolynomial{T,N}, c::S) where {T,N,S <: Number}
    R = eltype(one(T)/one(S))
    isinf(c)  && return zero(ImmutablePolynomial{R}, p.var)
    ImmutablePolynomial{R,N}(p.coeffs ./ c, p.var)
end

Base.:-(p::ImmutablePolynomial{T,N}) where {T,N} = ImmutablePolynomial{T,N}(.-p.coeffs, p.var)

Base.to_power_type(p::ImmutablePolynomial{T,N}) where {T,N} = p


## more performant versions borrowed from StaticArrays
## https://github.com/JuliaArrays/StaticArrays.jl/blob/master/src/linalg.jl
LinearAlgebra.norm(q::ImmutablePolynomial{T,0}) where {T} = zero(real(float(T)))
LinearAlgebra.norm(q::ImmutablePolynomial) = _norm(q.coeffs)
LinearAlgebra.norm(q::ImmutablePolynomial, p::Real) = _norm(q.coeffs, p)

@generated function _norm(a::NTuple{N,T}) where {T, N}

    expr = :(abs2(a[1]))
    for j = 2:N
        expr = :($expr + abs2(a[$j]))
    end
    
    return quote
        $(Expr(:meta, :inline))
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
        $(Expr(:meta, :inline))
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


