export   SparsePolynomial

"""
    SparsePolynomial{T, X}(coeffs::Dict, [var = :x])

Polynomials in the standard basis backed by a dictionary holding the
non-zero coefficients. For polynomials of high degree, this might be
advantageous. Addition and multiplication with constant polynomials
are treated as having no symbol.

# Examples:

```jldoctest
julia> using Polynomials

julia> P  = SparsePolynomial
SparsePolynomial

julia> p, q = P([1,2,3]), P([4,3,2,1])
(SparsePolynomial(1 + 2*x + 3*x^2), SparsePolynomial(4 + 3*x + 2*x^2 + x^3))

julia> p + q
SparsePolynomial(5 + 5*x + 5*x^2 + x^3)

julia> p * q
SparsePolynomial(4 + 11*x + 20*x^2 + 14*x^3 + 8*x^4 + 3*x^5)

julia> p + 1
SparsePolynomial(2 + 2*x + 3*x^2)

julia> q * 2
SparsePolynomial(8 + 6*x + 4*x^2 + 2*x^3)

julia> p = Polynomials.basis(P, 10^9) - Polynomials.basis(P,0) # also P(Dict(0=>-1, 10^9=>1))
SparsePolynomial(-1.0 + 1.0*x^1000000000)

julia> p(1)
0.0
```

"""
struct SparsePolynomial{T <: Number, X} <: StandardBasisPolynomial{T, X}
    coeffs::Dict{Int, T}
    function SparsePolynomial{T, X}(coeffs::AbstractDict{Int, S}) where {T <: Number, X, S}
        c = Dict{Int, T}(coeffs)
        for (k,v)  in coeffs
            iszero(v) && pop!(c,  k)
        end
        new{T, X}(c)
    end
    function SparsePolynomial{T,X}(checked::Val{false}, coeffs::AbstractDict{Int, T}) where {T <: Number, X}
        new{T,X}(coeffs)
    end
end

@register SparsePolynomial

function SparsePolynomial{T}(coeffs::AbstractDict{Int, S}, var::SymbolLike=:x) where {T <: Number, S}
    SparsePolynomial{T, Symbol(var)}(convert(Dict{Int,T}, coeffs))
end

function SparsePolynomial(coeffs::AbstractDict{Int, T}, var::SymbolLike=:x) where {T <: Number}
    SparsePolynomial{T, Symbol(var)}(coeffs)
end

function SparsePolynomial{T,X}(coeffs::AbstractVector{S}) where {T <: Number, X, S}

    if Base.has_offset_axes(coeffs)
        @warn "ignoring the axis offset of the coefficient vector"
    end

    offset = firstindex(coeffs)
    p = Dict{Int,T}(k - offset => v for (k,v) ∈ pairs(coeffs))
    return SparsePolynomial{T,X}(p)
end




# conversion
function Base.convert(P::Type{<:Polynomial}, q::SparsePolynomial)
    ⟒(P)(coeffs(q), indeterminate(q))
end

function Base.convert(P::Type{<:SparsePolynomial}, q::StandardBasisPolynomial{T}) where {T}
    R = promote_type(eltype(P), T)
    ⟒(P){R,indeterminate(P,q)}(coeffs(q))
end

## changes to common
degree(p::SparsePolynomial) = isempty(p.coeffs) ? -1 : maximum(keys(p.coeffs))
function isconstant(p::SparsePolynomial)
    n = length(keys(p.coeffs))
    (n > 1 || (n==1 && iszero(p[0]))) && return false
    return true
end

function basis(P::Type{<:SparsePolynomial}, n::Int)
    T,X = eltype(P), indeterminate(P)
    SparsePolynomial{T,X}(Dict(n=>one(T)))
end

# return coeffs as  a vector
# use p.coeffs to get Dictionary
function coeffs(p::SparsePolynomial{T})  where {T}

    n = degree(p)
    cs = zeros(T, n+1)
    for (k,v) in p.coeffs
        cs[k+1]=v
    end
    cs
    
end

# get/set index
function Base.getindex(p::SparsePolynomial{T}, idx::Int) where {T <: Number}
    get(p.coeffs, idx, zero(T))
end

function Base.setindex!(p::SparsePolynomial, value::Number, idx::Int)
    idx < 0  && return p
    if iszero(value)
        haskey(p.coeffs, idx) && pop!(p.coeffs, idx)
    else
        p.coeffs[idx]  = value
    end
    return p
end


Base.firstindex(p::SparsePolynomial) = sort(collect(keys(p.coeffs)), by=x->x[1])[1]
Base.lastindex(p::SparsePolynomial) = sort(collect(keys(p.coeffs)), by=x->x[1])[end]
Base.eachindex(p::SparsePolynomial) = sort(collect(keys(p.coeffs)), by=x->x[1])

# pairs iterates only over non-zero
# inherits order for underlying dictionary
function Base.iterate(v::PolynomialKeys{SparsePolynomial{T,X}}, state...) where {T,X}
    y = iterate(v.p.coeffs, state...)
    y == nothing && return nothing
    return (y[1][1], y[2])
end

function Base.iterate(v::PolynomialValues{SparsePolynomial{T,X}}, state...) where {T,X}
    y = iterate(v.p.coeffs, state...)
    y == nothing && return nothing
    return (y[1][2], y[2])
end



##
## ----
##
    
function evalpoly(x::S, p::SparsePolynomial{T}) where {T,S}

    tot = zero(T) * EvalPoly._one(x) 
    for (k,v) in p.coeffs
        tot = EvalPoly._muladd(x^k, v, tot)
    end
    
    return tot
    
end

#  map: over values -- not keys
function Base.map(fn, p::P, args...) where {P <: SparsePolynomial}
    ks, vs = keys(p.coeffs), values(p.coeffs)
    vs′ = map(fn, vs, args...)
    _convert(p, Dict(Pair.(ks, vs′)))
end


## Addition
function Base.:+(p::SparsePolynomial{T,X}, c::S) where {T, X, S <: Number}

    R = promote_type(T,S)
    P = SparsePolynomial{R,X}

    #D = Dict{Int, R}(kv for kv ∈ p.coeffs)
    D = Dict{Int, R}()
    for (k,v) ∈ pairs(p)
        @inbounds D[k] = v
    end
    @inbounds D[0] = get(D,0,zero(R)) + c
    iszero(D[0]) && pop!(D,0)

    return P(Val(false), D)
    
end

# Implement over fallback. A bit faster as it covers T != S
function Base.:+(p1::P1, p2::P2) where {T,X, P1<:SparsePolynomial{T,X},
                                        S,   P2<:SparsePolynomial{S,X}}

    R = promote_type(T,S)
    Q = SparsePolynomial{R,X}
    
    d1, d2 = degree(p1), degree(p2)
    cs = d1 > d2 ? ⊕(P1, p1.coeffs, p2.coeffs) : ⊕(P1, p2.coeffs, p1.coeffs)

    return d1 != d2 ? Q(Val(false), cs) : Q(cs)
    
end

## Multiplication
function Base.:*(p::P, c::S) where {T, X, P <: SparsePolynomial{T,X}, S <: Number}

    R = promote_type(T,S)
    Q = ⟒(P){R,X}
    
    q  = zero(Q)
    for (k,pₖ) ∈ pairs(p)
        q[k] = pₖ * c
    end
    
    return q
end

function Base.:*(p::P, q::Q) where {T,X,P<:SparsePolynomial{T,X},
                                    S,  Q<:SparsePolynomial{S,X}}
    R = promote_type(T,S)
    SparsePolynomial{R,X}(⊗(P, p.coeffs, q.coeffs))
end


function derivative(p::SparsePolynomial{T,X}, order::Integer = 1) where {T,X}
    
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return p

    R = eltype(one(T)*1)
    P = SparsePolynomial
    hasnan(p) && return P{R,X}(Dict(0 => R(NaN)))

    n = degree(p)
    order > n && return zero(P{R,X})

    dpn = zero(P{R,X})
    @inbounds for k in eachindex(p)
        dpn[k-order] =  reduce(*, (k - order + 1):k, init = p[k])
    end

    return dpn

end


function integrate(p::P) where {T, X, P<:SparsePolynomial{T,X}}
    
    R = eltype(one(T)/1)
    Q = SparsePolynomial{R,X}

    if hasnan(p)
        return Q(Dict(0 => NaN))
    end

    ∫p = zero(Q)
    for k in eachindex(p)
        ∫p[k + 1] = p[k] / (k+1)
    end

    return ∫p

end

