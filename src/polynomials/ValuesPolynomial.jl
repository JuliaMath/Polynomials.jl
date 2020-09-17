export ValuesPolynomial

"""
    ValuesPolynomial{T<:Number, N}(coeffs::AbstractVector{T}, var=:x)

Construct an immutable (static) polynomial from its coefficients `a`,
lowest order first, optionally in terms of the given variable `x`
where `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct
this through `ValuesPolynomial((a_0, a_1, ..., a_n))` (assuming
`a_n ≠ 0`). As well, a vector or number can be used for construction.

The usual arithmetic operators are overloaded to work with polynomials
as well as with combinations of polynomials and scalars. However,
operations involving two polynomials of different variables causes an
error, though for `+` and `*` operations, constant polynomials are
treated as having no variable. 

As the coefficient size is a compile-time constant, several performance
improvements are possible. For example, immutable polynomials can take advantage of 
faster polynomial evaluation provided by `evalpoly` from Julia 1.4.

    # Examples

```jldoctest
julia> using  Polynomials

julia> ValuesPolynomial((1, 0, 3, 4))
ValuesPolynomial(1 + 3*x^2 + 4*x^3)

julia> ValuesPolynomial((1, 2, 3), :s)
ValuesPolynomial(1 + 2*s + 3*s^2)

julia> one(ValuesPolynomial)
ValuesPolynomial(1.0)
```

!!! note
    This was modeled after https://github.com/tkoolen/StaticUnivariatePolynomials.jl by @tkoolen.

"""
struct ValuesPolynomial{T <: Number,  N} <: StandardBasisPolynomial{T}
    coeffs::Values{N, T}
    var::Symbol
    
    function ValuesPolynomial{T,N}(coeffs::Values{N,T}, var::SymbolLike=:x) where {T <: Number,N}
        N == 0 && return new{T,0}(coeffs, Symbol(var))
        iszero(coeffs[end]) &&  throw(ArgumentError("Leading  term must  be  non-zero"))
        new{T,N}(coeffs, Symbol(var))
    end

    function ValuesPolynomial{T,N}(coeffs::NTuple{N,T}, var::SymbolLike=:x) where {T <: Number, N}
        new{T, N}(Values{N,T}(coeffs), var)
    end

    function ValuesPolynomial{T,N}(coeffs::AbstractVector{S}, var::SymbolLike=:x) where {T <: Number, N, S}
        new{T,N}(Values{N,T}(tuple(coeffs...)), var)
    end

    function ValuesPolynomial{T}(coeffs::Values{M,S}, var::SymbolLike=:x) where {T, S<: Number, M}
        N = findlast(!iszero, coeffs)
        if N == nothing
            return zero(ValuesPolynomial{T}, var)
        else
            cs = NTuple{N,T}(coeffs[i] for  i in  1:N)
        end
        new{T,N}(Values(cs), var)
    end
    
    function ValuesPolynomial{T}(coeffs::NTuple{M,S}, var::SymbolLike=:x) where {T, S<: Number, M}
        ValuesPolynomial{T}(Values(coeffs), var)
    end

    # entry point from abstract.jl; note Vector type
    function ValuesPolynomial{T}(coeffs::Vector{T}, var::SymbolLike=:x) where {T}
        M = length(coeffs)
        ValuesPolynomial{T}(Values{M,T}(tuple(coeffs...)), var)
    end
    
end

@register ValuesPolynomial

# less specific than NTuple
function ValuesPolynomial{T,N}(coeffs::Tuple, var::SymbolLike=:x)  where {T,N}
    ValuesPolynomial{T,N}(T.(coeffs), var)
end

function ValuesPolynomial{T}(coeffs::Tuple, var::SymbolLike=:x)  where {T}
    ValuesPolynomial{T}(T.(coeffs), var)
end

function ValuesPolynomial(coeffs::Tuple, var::SymbolLike=:x)
    cs = NTuple(promote(coeffs...))
    T = eltype(cs)
    ValuesPolynomial{T}(cs, var)
end


# function ValuesPolynomial(coeffs::Values{M,T}, var::SymbolLike=:x) where {M, T}
#     ValuesPolynomial{T}(coeffs, var)
# end



# Convenience; pass tuple to Polynomial
# Not documented, not sure this is a good idea as P(...)::P is not true...
Polynomial(coeffs::Values{N,T}, var::SymbolLike = :x) where{N,T} =
    ValuesPolynomial(coeffs, var)

function Polynomial{T}(coeffs::Values{N,S}, var::SymbolLike = :x) where{N,T,S}
    ValuesPolynomial{N,T}(T.(coeffs), var)
end

##
## ----
##
# overrides from common.jl due to  coeffs being non mutable, ...

Base.collect(p::P) where {P <: ValuesPolynomial} = [pᵢ for pᵢ ∈ p]

Base.copy(p::P) where {P <: ValuesPolynomial} = P(coeffs(p), p.var)

# catch q == 0 case

#  zero, one, variable
function Base.zero(P::Type{<:ValuesPolynomial}, var::SymbolLike=:x)
    R = eltype(P)
    ValuesPolynomial{R,0}(Values{0,R}(), var)
end

function  Base.one(P::Type{<:ValuesPolynomial}, var::SymbolLike=:x)
    R = eltype(P)
    ValuesPolynomial{R,1}(Values{1,R}(1),var)
end
function variable(P::Type{<:ValuesPolynomial}, var::SymbolLike=:x)
    R  = eltype(P)
    ValuesPolynomial{R,2}(Values{2,R}((0,1)),var)
end

# degree, isconstant
degree(p::ValuesPolynomial{T,N}) where {T,N} = N - 1
isconstant(p::ValuesPolynomial{T,N}) where {T,N}  = N <= 1

function Base.getindex(p::ValuesPolynomial{T,N}, idx::Int) where {T <: Number,N}
    (idx <  0 || idx > N-1) && return zero(T)
    return p.coeffs[idx + 1]
end

Base.setindex!(p::ValuesPolynomial, val::Number,  idx::Int) = throw(ArgumentError("Values polynomials are immutable"))

for op in [:isequal, :(==)]
    @eval function Base.$op(p1::ValuesPolynomial{T,N}, p2::ValuesPolynomial{S,M}) where {T,N,S,M}
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
## This fixes the last two cases for P=ValuesPolynomial, and could replace
## isapprox in common.jl
##
## check Inf values (matching Vector ≈ Vector)
## @test P([Inf]) ≈ Inf
## @test !(P([Inf]) ≈ P([-Inf])) # default compares zero(P) to zero(P)
## @test !(P([1,2,Inf]) ≈ P([1,3,Inf])) # default uses truncate with  Inf which clobbers numbers
function Base.isapprox(p1::ValuesPolynomial{T},
    p2::ValuesPolynomial{S};
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
                       atol::Real = 0,) where {T,S, P<:ValuesPolynomial{T}}
    return isapprox(p1, ⟒(P){T}(n,p1.var))
end

function Base.chop(p::ValuesPolynomial{T,N};
              rtol::Real = Base.rtoldefault(real(T)),
              atol::Real = 0)  where {T,N}
    cs = coeffs(p)
    for i in N:-1:1
        if !isapprox(cs[i], zero(T), rtol=rtol, atol=atol)
            return ValuesPolynomial{T,i}(cs[1:i], p.var)
        end
    end
    zero(ValuesPolynomial{T}, p.var)
end

function Base.truncate(p::ValuesPolynomial{T,N};
                       rtol::Real = Base.rtoldefault(real(T)),
                       atol::Real = 0)  where {T,N}
    q = chop(p, rtol=rtol, atol=atol)
    iszero(q) && return  q
    cs = coeffs(q)
    thresh = maximum(abs,cs) * rtol + atol
    cs′ = map(c->abs(c) <= thresh ? zero(T) : c, cs)
    ValuesPolynomial{T}(tuple(cs′...), p.var)
end

# no in-place chop! and truncate!
chop!(p::ValuesPolynomial; kwargs...) =  chop(p; kwargs...)
truncate!(p::ValuesPolynomial; kwargs...) =  truncate(p; kwargs...)

LinearAlgebra.conj(p::P) where {P <: ValuesPolynomial} = P(conj([aᵢ for aᵢ in coeffs(p)]))


##
## --------------------
##

(p::ValuesPolynomial{T,N})(x::S) where {T,N,S} = evalpoly(x, p.coeffs.v)


function Base.:+(p1::ValuesPolynomial{T,N}, p2::ValuesPolynomial{S,M}) where {T,N,S,M}

    R = promote_type(S,T)
    iszero(N) && return ValuesPolynomial{R}(coeffs(p2), p2.var)
    iszero(M) && return ValuesPolynomial{R}(coeffs(p1), p1.var)
    
    isconstant(p1) && p1.var != p2.var && return p2 + p1[0]*one(ValuesPolynomial{T},p2.var)
    isconstant(p2) && p1.var != p2.var && return p1 + p2[0]*one(ValuesPolynomial{S},p1.var)

    p1.var != p2.var && error("Polynomials must have same variable")


    if  N == M
        cs = NTuple{N,R}(p1[i] + p2[i] for i in 0:N-1)
        ValuesPolynomial{R}(cs, p1.var)        
    elseif N < M
        cs = (p2.coeffs.v) ⊕ (p1.coeffs.v)
        ValuesPolynomial{R,M}(cs, p1.var)        
    else
        cs = (p1.coeffs.v) ⊕ (p2.coeffs.v)
        ValuesPolynomial{R,N}(cs, p1.var)                
    end

    

end


function Base.:*(p1::ValuesPolynomial{T,N}, p2::ValuesPolynomial{S,M}) where {T,N,S,M}
    isconstant(p1) && return p2 * p1[0] 
    isconstant(p2) && return p1 * p2[0]
    p1.var != p2.var && error("Polynomials must have same variable")
    R = promote_type(S,T)
    cs = (p1.coeffs.v) ⊗ (p2.coeffs.v)
    ValuesPolynomial{R, N+M-1}(cs, p1.var)
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
function Base.:+(p::ValuesPolynomial{T,N}, c::S) where {T, N, S<:Number}
    R = promote_type(T,S)

    iszero(c) && return ValuesPolynomial{R,N}(p.coeffs, p.var)
    N == 0 && return ValuesPolynomial{R,1}((c,), p.var)
    N == 1 && return ValuesPolynomial((p[0]+c,), p.var)

    q = ValuesPolynomial{R,1}((c,), p.var)
    return p + q

end

function Base.:*(p::ValuesPolynomial{T,N}, c::S) where {T, N, S <: Number}
    R = promote_type(T,S)
    iszero(c) && return zero(ValuesPolynomial{R}, p.var)
    ValuesPolynomial{R,N}(p.coeffs .* c, p.var)
end

function Base.:/(p::ValuesPolynomial{T,N}, c::S) where {T,N,S <: Number}
    R = eltype(one(T)/one(S))
    isinf(c)  && return zero(ValuesPolynomial{R}, p.var)
    ValuesPolynomial{R,N}(p.coeffs ./ c, p.var)
end

Base.:-(p::ValuesPolynomial{T,N}) where {T,N} = ValuesPolynomial{T,N}(-p.coeffs, p.var)

Base.to_power_type(p::ValuesPolynomial{T,N}) where {T,N} = p
