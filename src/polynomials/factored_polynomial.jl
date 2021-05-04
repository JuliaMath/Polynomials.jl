export FactoredPolynomial

"""
    FactoredPolynomial{T,X}

A polynomial type that stores its data in a dictionary whose keys are the roots and whose values are the respecitive multiplicities along with a leading coefficient.

The structure is utilized for polynomial multiplication and powers, the finding of roots, and the identification of a greatest common divisor.

## Examples

```jldoctest factored_polynomial
julia> using Polynomials

julia> p = FactoredPolynomial(Dict([0=>1, 1=>2, 3=>4]))
FactoredPolynomial(1 * x * (x - 3)⁴ * (x - 1)²)

julia> q = fromroots(FactoredPolynomial, [0,1,2,3])
FactoredPolynomial(1 * x * (x - 2) * (x - 3) * (x - 1))

julia> p*q
FactoredPolynomial(1 * x² * (x - 2) * (x - 3)⁵ * (x - 1)³)

julia> p^1000
FactoredPolynomial(1 * x¹⁰⁰⁰ * (x - 3)⁴⁰⁰⁰ * (x - 1)²⁰⁰⁰)

julia> gcd(p,q)
FactoredPolynomial(1 * x * (x - 3) * (x - 1))
```
"""
struct FactoredPolynomial{T <: Number, X} <: StandardBasisPolynomial{T, X}
    coeffs::Dict{T,Int}
    c::T
    function FactoredPolynomial{T, X}(cs::Dict{T,Int}, c=one(T)) where {T, X}
        D = Dict{T,Int}()
        for (k,v) ∈ cs
            v > 0 && (D[k] = v)
        end
        new{T,X}(D,T(c))
    end
    function FactoredPolynomial(cs::Dict{T,Int}, c::S=1, var::SymbolLike=:x) where {T,S}
        X = Symbol(var)
        R = promote_type(T,S)
        D = convert(Dict{R,Int}, cs)
        FactoredPolynomial{R,X}(D, R(c))
    end
end

function FactoredPolynomial{T,X}(coeffs::AbstractVector{S}) where {T,S,X}
    p = Polynomial{T,X}(Val(false), coeffs)
    zs = Multroot.multroot(p)
    c = p[end]
    D = Dict(zip(zs.values, zs.multiplicities))
    FactoredPolynomial(D, c)
end

function FactoredPolynomial(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {T}
    X = Symbol(var)
    p = Polynomial{T,X}(Val(false), coeffs)
    zs = Multroot.multroot(p)
    c = p[end]
    D = Dict(zip(zs.values, zs.multiplicities))
    FactoredPolynomial(D, c)
end


## ----

function printpoly(io::IO, p::FactoredPolynomial{T,X}, mimetype=nothing) where {T,X}
    if iszero(p.c)
        print(io, p.c)
    else
        isone(p.c) && iszero(length(p.coeffs)) && print(io, p.c)
        !isone(p.c) && print(io, p.c)
        x = string(X)
        for (i,(k,v)) ∈ enumerate(p.coeffs)
            (!isone(p.c) || i != 1) && print(io, " * ")
            if iszero(k)
                print(io, x)
            else
                print(io, "(")
                print(io, x)
                if hasneg(T)
                    isneg(k) ?  print(io, " + ", -k) : print(io, " - ", k) 
                else
                    print(io, " - ", k)
                end
                print(io, ")")
            end
            v > 1 && unicode_exponent(io, v)
        end
    end
end

Base.zero(::Type{FactoredPolynomial{T,X}}) where {T, X} = FactoredPolynomial{T,X}(Dict{T,Int}(), one(T))
Base.one(::Type{FactoredPolynomial{T,X}}) where {T, X} = FactoredPolynomial{T,X}(Dict{T,Int}(), one(T))
variable(::Type{FactoredPolynomial{T,X}}) where {T, X} = FactoredPolynomial{T,X}(Dict{T,Int}(0=>1), one(T))


# abstract
Base.convert(::Type{P}, p::P) where {P <: FactoredPolynomial} = p
function Base.convert(P::Type{<:FactoredPolynomial}, p::FactoredPolynomial{T,X}) where {T,X}
    T′ = _eltype(P)
    𝑻 = T′ == nothing ? T : T′
    d = Dict{𝑻,Int}()
    copy!(d, p.coeffs)
    FactoredPolynomial{𝑻,X}(d, p.c)
end
Base.promote(p::P,q::Q) where {X,T,P<:FactoredPolynomial{T,X},Q<:FactoredPolynomial{T,X}} = p,q
Base.promote_rule(::Type{<:FactoredPolynomial{T,X}}, ::Type{<:FactoredPolynomial{S,X}}) where {T,S,X} =
    FactoredPolynomial{promote_type(T,S), X}
Base.promote_type(::Type{<:FactoredPolynomial{T,X}}, ::Type{<:FactoredPolynomial{S,X}}) where {T,S,X} =
    FactoredPolynomial{promote_type(T,S), X}
Base.promote_type(::Type{<:FactoredPolynomial{T,X}}, ::Type{S}) where {T,S,X} =
    FactoredPolynomial{promote_type(T,S), X}
FactoredPolynomial{T,X}(n::S) where {T,X,S<:Number} = T(n) * one(FactoredPolynomial{T,X})
FactoredPolynomial{T}(n::S, var::SymbolLike=:x) where {T,S} = T(n) * one(FactoredPolynomial{T,X})
FactoredPolynomial(n::S, var::SymbolLike=:x) where {S<:Number} = n * one(FactoredPolynomial{S,Symbol(var)})
FactoredPolynomial(var::SymbolLike=:x) = variable(FactoredPolynomial, Symbol(var))
(p::FactoredPolynomial)(x) = evalpoly(x, p)

function Base.convert(::Type{Polynomial}, p::FactoredPolynomial{T,X}) where {T,X}
    x = variable(Polynomial{T,X})
    isconstant(p) && return Polynomial{T,X}(p.c)
    p(x)
end
function Base.convert(P::Type{FactoredPolynomial}, p::Polynomial{T,X}) where {T,X}
    isconstant(p) && return P(constantterm(p), X)
    P(coeffs(p), X)
end
    
#
function evalpoly(x, p::FactoredPolynomial)
    iszero(length(p.coeffs)) && return p.c
    p.c * prod((x-k)^v for (k,v) ∈ p.coeffs)
end


function fromroots(::Type{P}, r::AbstractVector{T}; var::SymbolLike=:x) where {T <: Number, P<:FactoredPolynomial}
    d = Dict{T,Int}()
    for rᵢ ∈ r
        d[rᵢ] = get(d, rᵢ, 0) + 1
    end
    FactoredPolynomial{T, Symbol(var)}(d)
end

roots(p::FactoredPolynomial) = collect(keys(p.coeffs))

coeffs(p::FactoredPolynomial) = coeffs(convert(Polynomial, p))

function degree(p::FactoredPolynomial)
    d = length(p.coeffs)
    d > 0 && return sum(values(p.coeffs))
    iszero(p.c)  ? -1 : 0
end

function integrate(p::FactoredPolynomial)
    error("XXX")
end

function derivative(p::FactoredPolynomial,n::Int)
    error("XXX")
end

# pairs,keys, values
Base.keys(p::FactoredPolynomial) = keys(convert(Polynomial, p))
Base.values(p::FactoredPolynomial) = values(convert(Polynomial, p))
Base.pairs(p::FactoredPolynomial) = pairs(convert(Polynomial, p))


# addition 
function Base.:+(p::P, q::P) where {T,X,P<:FactoredPolynomial{T,X}}
    convert(P, convert(Polynomial, p) + convert(Polynomial, q))
end

# multiplication
function Base.:*(p::P, q::P) where {T,X, P<:FactoredPolynomial{T,X}}
    d = copy(p.coeffs)
    for (k,v) ∈ q.coeffs
        if haskey(d,k)
            d[k] += v
        else
            d[k] = v
        end
    end
    P(d, p.c*q.c)
end

function Base.:*(p::P, c::S) where {S, T, X, P <: FactoredPolynomial{T, X}}
    R = promote_type(T,S)
    d = Dict{R,Int}()
    for (k,v) ∈ p.coeffs
        d[k] = v
    end
    FactoredPolynomial{R,X}(d, c * p.c)
end

function Base.:/(p::P, c::S) where {S, T, X, P <: FactoredPolynomial{T, X}}
    p * (1/c)
end

function Base.:^(p::FactoredPolynomial{T,X}, n::Integer) where {T,X}
    d = Dict{T,Int}()
    for (k,v) ∈ p.coeffs
        d[k] = v*n
    end
    FactoredPolynomial{T,X}(d, p.c^n)
end

## gcd, divrem, uvw
function Base.gcd(p::FactoredPolynomial{T,X}, q::FactoredPolynomial{T,X}) where {T, X}
    d = Dict{T,Int}()

    for k ∈ intersect(keys(p.coeffs), keys(q.coeffs))
        d[k] = min(p.coeffs[k], q.coeffs[k])
    end
    FactoredPolynomial{T,X}(d)
end

# return u,v,w with p = u*v , q = u*w
function uvw(p::FactoredPolynomial{T,X}, q::FactoredPolynomial{T,X}) where {T, X}
    du, dv, dw = Dict{T,Int}(), Dict{T,Int}(), Dict{T,Int}()
    dp,dq = p.coeffs, q.coeffs
    kp,kq = keys(dp), keys(dq)
    
    for k ∈ setdiff(kp, kq)
        dv[k] = dp[k]
    end
    for k ∈ setdiff(kq, kp)
        dw[k] = dq[k]
    end
    for k ∈ intersect(kp, kq)
        pₖ,qₖ = dp[k], dq[k]
        m = min(pₖ, qₖ)
        du[k] = m
        dv[k] = pₖ - m; 
        dw[k] = qₖ - m
    end
    FactoredPolynomial(du), FactoredPolynomial(dv, p.c), FactoredPolynomial(dw, q.c)
end
        
        

function Base.divrem(p::FactoredPolynomial{T,X}, q::FactoredPolynomial{T,X}) where {T, X}
    u,v,w = uvw(p,q)
    isconstant(w) && return (v / q.c, zero(v))
    vv, ww = convert(Polynomial, v), convert(Polynomial, w)
    d,r = divrem(vv,ww)
    dd, rr = convert(FactoredPolynomial, d), convert(FactoredPolynomial,r)
    dd,rr
end    
