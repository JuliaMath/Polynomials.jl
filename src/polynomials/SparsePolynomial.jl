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
    function SparsePolynomial{T, X}(coeffs::AbstractDict{Int, T}) where {T <: Number, X}
        c = Dict(coeffs)
        for (k,v)  in coeffs
            iszero(v) && pop!(c,  k)
        end
        new{T, X}(c)
    end
end

@register SparsePolynomial

function SparsePolynomial{T,X}(coeffs::AbstractVector{T}) where {T <: Number, X}
    firstindex(coeffs) >= 0 || throw(ArgumentError("Use the `LaurentPolynomial` type for arrays with a negative first index"))

    if Base.has_offset_axes(coeffs)
      @warn "ignoring the axis offset of the coefficient vector"
    end
    c = OffsetArrays.no_offset_view(coeffs) # ensure 1-based indexing
    p = Dict{Int,T}(i - 1 => v for (i,v) in pairs(c))
    return SparsePolynomial{T,X}(p)
end

function SparsePolynomial(coeffs::AbstractDict{Int, T}, var::SymbolLike=:x) where {T <: Number}
    SparsePolynomial{T, Symbol(var)}(coeffs)
end



# conversion
function Base.convert(P::Type{<:Polynomial}, q::SparsePolynomial)
    ⟒(P)(coeffs(q), var(q))
end

function Base.convert(P::Type{<:SparsePolynomial}, q::StandardBasisPolynomial{T}) where {T}
    R = promote(eltype(P), T)
    ⟒(P){R}(coeffs(q), var(q))
end

## changes to common
degree(p::SparsePolynomial) = isempty(p.coeffs) ? -1 : maximum(keys(p.coeffs))
function isconstant(p::SparsePolynomial)
    n = length(keys(p.coeffs))
    (n > 1 || (n==1 && iszero(p[0]))) && return false
    return true
end

function basis(P::Type{<:SparsePolynomial}, n::Int, var::SymbolLike=:x) 
    T = eltype(P)
    X = Symbol(var)
    SparsePolynomial{T,X}(Dict(n=>one(T)))
end

# return coeffs as  a vector
# use p.coeffs to get Dictionary
function  coeffs(p::SparsePolynomial{T})  where {T}

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

# only from tail
function chop!(p::SparsePolynomial{T};
               rtol::Real = Base.rtoldefault(real(T)),
               atol::Real = 0,) where {T}

    for k in sort(collect(keys(p.coeffs)), by=x->x[1], rev=true)
        if isapprox(p[k], zero(T); rtol = rtol, atol = atol)
            pop!(p.coeffs, k)
        else
            return p
        end
    end
    
    return p
    
end

function truncate!(p::SparsePolynomial{T};
                   rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {T}
    
    max_coeff = maximum(abs, coeffs(p))
    thresh = max_coeff * rtol + atol

    for (k,val) in  p.coeffs
        if abs(val) <= thresh
            pop!(p.coeffs,k)
        end
    end
    
    return p
    
end

##
## ----
##
    

function (p::SparsePolynomial{T})(x::S) where {T,S}
    
    tot = zero(T) * _one(x) 
    for (k,v) in p.coeffs
        tot = _muladd(x^k, v, tot)
    end
    
    return tot
    
end

#  map: over values -- not keys
function Base.map(fn, p::P, args...) where {P <: SparsePolynomial}
    ks, vs = keys(p.coeffs), values(p.coeffs)
    vs′ = map(fn, vs, args...)
    _convert(p, Dict(Pair.(ks, vs′)))
end


   
function Base.:+(p1::SparsePolynomial{T,X}, p2::SparsePolynomial{S,Y}) where {T, X, S, Y}

    isconstant(p1) && return p2 + p1[0]
    isconstant(p2) && return p1 + p2[0]

    X != Y && error("SparsePolynomials must have same variable")

    R = promote_type(T,S)
    p = zero(SparsePolynomial{R,X})

    # this allocates in the union
#    for i in union(eachindex(p1), eachindex(p2)) 
#        p[i] = p1[i] + p2[i]
#    end

    # this seems faster
    for i in eachindex(p1)
        p[i] = p1[i] + p2[i]
    end
    for i in eachindex(p2)
        if iszero(p[i])
            @inbounds p[i] = p1[i] + p2[i]
        end
    end
    

    return  p

end

function Base.:+(p::SparsePolynomial{T,X}, c::S) where {T, X, S <: Number}

    R = promote_type(T,S)
    P = SparsePolynomial
    
    q = zero(P{R,X})
    for k in eachindex(p)
        @inbounds q[k] = R(p[k])
    end
    q[0] = q[0] + c

    return q
end

function Base.:*(p1::SparsePolynomial{T,X}, p2::SparsePolynomial{S,Y}) where {T,X,S,Y}

    isconstant(p1) && return p2 * p1[0]
    isconstant(p2) && return p1 * p2[0]
    X != Y && error("SparsePolynomials must have same variable")

    R = promote_type(T,S)
    P = SparsePolynomial
    
    p  = zero(P{R, X})
    for i in eachindex(p1)
        p1ᵢ = p1[i]
        for j in eachindex(p2)
            @inbounds p[i+j] = muladd(p1ᵢ, p2[j], p[i+j])
        end
    end
    
    return p
    
end


function Base.:*(p::P, c::S) where {T, X, P <: SparsePolynomial{T,X}, S <: Number}

    R = promote_type(T,S)
    q  = zero(⟒(P){R,X})
    for k in eachindex(p)
        q[k] = p[k] * c
    end
    
    return q
end



function derivative(p::SparsePolynomial{T,X}, order::Integer = 1) where {T,X}
    
    order < 0 && error("Order of derivative must be non-negative")
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


function integrate(p::SparsePolynomial{T,X}, k::S) where {T, X, S<:Number}
    
    R = eltype((one(T)+one(S))/1)
    P = SparsePolynomial

    if hasnan(p) || isnan(k)
        return P{R,X}(Dict(0 => R(NaN))) # not R(NaN)!! don't like XXX
    end

    ∫p = P{R,X}(R(k))
    for k in eachindex(p)
        ∫p[k + 1] = p[k] / (k+1)
    end
    
    return ∫p
    
end
