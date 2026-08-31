export FactoredPolynomial

"""
    FactoredPolynomial{T,X}

A polynomial type that stores its data in a dictionary whose keys are the roots and whose values are the respective multiplicities along with a leading coefficient.

The structure is utilized for scalar multiplication, polynomial multiplication and powers, the finding of roots, and the identification of a greatest common divisor. For other operations, say addition, the operation is done after converting to the `Polynomial` type then converting back. (This requires the identification of the roots, so is subject to numeric issues.)

## Examples

```jldoctest factored_polynomial
julia> using Polynomials

julia> p = FactoredPolynomial(Dict([0=>1, 1=>2, 3=>4]))
FactoredPolynomial(x * (x - 1)² * (x - 3)⁴)

julia> q = fromroots(FactoredPolynomial, [0,1,2,3])
FactoredPolynomial(x * (x - 1) * (x - 2) * (x - 3))

julia> p*q
FactoredPolynomial(x² * (x - 1)³ * (x - 2) * (x - 3)⁵)

julia> p^1000
FactoredPolynomial(x¹⁰⁰⁰ * (x - 1)²⁰⁰⁰ * (x - 3)⁴⁰⁰⁰)

julia> gcd(p,q)
FactoredPolynomial(x * (x - 1) * (x - 3))

julia> p = Polynomial([24, -50, 35, -10, 1])
Polynomial(24 - 50*x + 35*x^2 - 10*x^3 + x^4)

julia> q = convert(FactoredPolynomial, p); # noisy form of `factor`, subject to floating point issues

julia> map(x->round(x, digits=12), q) # map works over factors and leading coefficient -- not coefficients in the standard basis
FactoredPolynomial((x - 1.0) * (x - 2.0) * (x - 3.0) * (x - 4.0))
```
"""
struct FactoredPolynomial{T <: Number, X} <: AbstractPolynomial{T, X}
    coeffs::OrderedDict{T,Int}
    c::T
    function FactoredPolynomial{T, X}(checked::Val{false}, cs::AbstractDict{T,Int}, c::T) where {T, X}
        #new{T,X}(convert(OrderedDict,cs),T(c))
        D = OrderedDict(sort(collect(pairs(cs)); by = x -> x.first, rev=false))
        new{T,X}(D,T(c))
    end
    function FactoredPolynomial{T, X}(cs::AbstractDict{T,Int}, c=one(T)) where {T, X}
        D = OrderedDict{T,Int}()
        for (k,v) ∈ cs
            v > 0 && (D[k] = v)
        end
        FactoredPolynomial{T,X}(Val(false), D,T(c))
    end
    function FactoredPolynomial(cs::AbstractDict{T,Int}, c::S=1, var::SymbolLike=:x) where {T,S}
        X = Symbol(var)
        R = promote_type(T,S)
        D = OrderedDict{R,Int}(sort(collect(pairs(cs)); by = x -> x.first, rev=false))
        FactoredPolynomial{R,X}(D, R(c))
    end
end

# There are idiosyncrasies with this type
# * unlike P{T}(...) we allow T to widen here if the roots of the polynomial Polynomial(coeffs) needs
#   a wider type (e.g. Complex{Float64}, not Float64)
# * the handling of Inf and NaN, when a specified coefficient, is to just return a constant poly (Inf or NaN)
function FactoredPolynomial{T,X}(coeffs::AbstractVector{S}) where {T,S,X}
    p = Polynomial{T,X}(T.(coeffs))
    P = FactoredPolynomial{T,X}
    iszero(p)     && return zero(P)
    isconstant(p) && return FactoredPolynomial{T,X}(Val(false), Dict{T,Int}(), T(coeffs[1]))
    any(isnan, p) && return FactoredPolynomial{T,X}(Val(false), Dict{T,Int}(), T(NaN))
    any(isinf, p) && return FactoredPolynomial{T,X}(Val(false), Dict{T,Int}(), T(Inf))

    zs = Multroot.multroot(p)
    c = p[end]

    D = LittleDict(k=>v for (k,v) ∈ zip(zs.values, zs.multiplicities))
    D′ = OrderedDict{eltype(zs.values), Int}()
    for z ∈ lejaorder(zs.values)
        D′[z] = D[z]
    end

    FactoredPolynomial(D′, c, X)
end

function FactoredPolynomial{T}(coeffs::AbstractVector{S}, var::SymbolLike=:x) where {T,S}
    X = Symbol(var)
    FactoredPolynomial{T,X}(coeffs)
end

function FactoredPolynomial(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {T}
    X = Symbol(var)
    FactoredPolynomial{T,X}(coeffs)
end



## ----
# abstract.jl The use of @register FactoredPolynomial didn't quite work, so
# that is replicated here and modified.
Base.convert(::Type{P}, p::P) where {P <: FactoredPolynomial} = p
function Base.convert(P::Type{<:FactoredPolynomial}, p::FactoredPolynomial{T,X}) where {T,X}
    T′ = _eltype(P)
    𝑻 = isnothing(T′) ? T : T′
    𝑿 = indeterminate(P, p)
    d = Dict{𝑻,Int}()
    copy!(d, p.coeffs)
    FactoredPolynomial{𝑻,𝑿}(d, p.c)
end

Base.promote(p::P,q::Q) where {X,T,P<:FactoredPolynomial{T,X},Q<:FactoredPolynomial{T,X}} = p,q

Base.promote_rule(::Type{<:FactoredPolynomial{T,X}}, ::Type{<:FactoredPolynomial{S,X}}) where {T,S,X} =
    FactoredPolynomial{promote_type(T,S), X}

Base.promote_rule(::Type{<:FactoredPolynomial{T,X}}, ::Type{S}) where {T,S<:Number,X} =
    FactoredPolynomial{promote_type(T,S), X}

FactoredPolynomial{T,X}(n::S) where {T,X,S<:Number} = T(n) * one(FactoredPolynomial{T,X})

FactoredPolynomial{T}(n::S, var::SymbolLike=:x) where {T,S<:Number} = T(n) * one(FactoredPolynomial{T,Symbol(var)})

FactoredPolynomial(n::S, var::SymbolLike=:x) where {S<:Number} = n * one(FactoredPolynomial{S,Symbol(var)})

FactoredPolynomial(var::SymbolLike=:x) = variable(FactoredPolynomial, Symbol(var))

(p::FactoredPolynomial)(x) = evalpoly(x, p)

function Base.convert(::Type{<:Polynomial}, p::FactoredPolynomial{T,X}) where {T,X}
    x = variable(Polynomial{T,X})
    isconstant(p) && return Polynomial{T,X}(p.c)
    p(x)
end

function Base.convert(P::Type{<:FactoredPolynomial}, p::Polynomial{T,X}) where {T,X}
    isconstant(p) && return ⟒(P)(constantterm(p), X)
    ⟒(P)(coeffs(p), X)
end

function copy_with_eltype(::Type{T}, ::Val{X}, p::P) where {T, X, S, Y, P<:FactoredPolynomial{S, Y}}
    d = convert(Dict{T, Int}, p.coeffs)
    FactoredPolynomial{T, X}(d)
end

## ----
## apply map to factors and the leading coefficient, not the coefficients
function Base.map(fn, p::P, args...)  where {T,X,P<:FactoredPolynomial{T,X}}
    𝒅 = Dict{T, Int}()
    for (k,v) ∈ p.coeffs
        𝒌 = fn(k, args...)
        𝒅[𝒌] = v
    end
    𝒄 = fn(p.c, args...)
    P(𝒅,𝒄)
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

## ----
Base.lastindex(p::FactoredPolynomial) = degree(p)
Base.eachindex(p::FactoredPolynomial) = 0:degree(p)

function Base.getindex(p::FactoredPolynomial{T}, idx::Int) where {T <: Number}
    m,M = firstindex(p), lastindex(p)
    idx < m && throw(BoundsError(p, idx))
    idx > M && return zero(T)
    coeffs(p)[idx - m + 1]
end

Base.setindex!(p::FactoredPolynomial, val,  idx::Int) = throw(ArgumentError("FactoredPolynomials are immutable"))

# pairs,keys, values
Base.keys(p::FactoredPolynomial)   = keys(convert(Polynomial, p))
Base.values(p::FactoredPolynomial) = values(convert(Polynomial, p))
Base.pairs(p::FactoredPolynomial)  = pairs(convert(Polynomial, p))


Base.copy(p::P) where {P<:FactoredPolynomial} = P(copy(p.coeffs), p.c)
Base.:(==)(p1::FactoredPolynomial, p2::FactoredPolynomial) =
    check_same_variable(p1,p2) && p1.c == p2.c && (p1.coeffs == p2.coeffs)

# what does it mean to be approximate here?
# Question pushed off to the Polynomial type
# might be better to compare roots and multiplicities
function Base.isapprox(p1::FactoredPolynomial{T,X},
                       p2::FactoredPolynomial{S,Y};
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {T,X,S,Y}
    assert_same_variable(p1, p2)

    # p ≈ q if 𝒑 ≈ 𝒒
    𝑷 = Polynomial
    𝒑,𝒒 = convert(𝑷,p1), convert(𝑷,p2)
    return isapprox(𝒑, 𝒒, atol=atol, rtol=rtol)

end


## ----
domain(::Type{<:FactoredPolynomial}) = Interval{Open,Open}(-Inf, Inf)
mapdomain(::Type{<:FactoredPolynomial}, x::AbstractArray) = x
Base.iszero(p::FactoredPolynomial) = iscoeffzero(p.c)

Base.zero(::Type{FactoredPolynomial{T,X}}) where {T, X} = FactoredPolynomial{T,X}(Dict{T,Int}(),     zero(T))
Base.one(::Type{FactoredPolynomial{T,X}}) where {T, X}  = FactoredPolynomial{T,X}(Dict{T,Int}(),     one(T))
variable(::Type{FactoredPolynomial{T,X}}) where {T, X}  = FactoredPolynomial{T,X}(Dict{T,Int}(0=>1), one(T))


#
function evalpoly(x, p::FactoredPolynomial)
    iszero(length(p.coeffs)) && return p.c * one(x)
    p.c * prod((x - k)^v for (k,v) ∈ p.coeffs)
end

function evalpoly(x::Array, p::FactoredPolynomial)
    iszero(length(p.coeffs)) && return p.c * one(x)
    p.c * prod((x - k*I)^v for (k,v) ∈ p.coeffs)
end

coeffs(p::FactoredPolynomial) = coeffs(convert(Polynomial, p))

function degree(p::P) where {P <: FactoredPolynomial}
    d = length(p.coeffs)
    d > 0 && return sum(values(p.coeffs))
    iszero(p.c)  ? -1 : 0
end

## ----

function fromroots(::Type{P}, r::AbstractVector{T}; var::SymbolLike=:x) where {T <: Number, P<:FactoredPolynomial}
    X = Symbol(var)
    d = OrderedDict{T,Int}()
    for rᵢ ∈ lejaorder(r)
        d[rᵢ] = get(d, rᵢ, 0) + 1
    end
    FactoredPolynomial{T, X}(d)
end

roots(p::FactoredPolynomial{T}) where {T} = Base.typed_vcat(T,[repeat([k],v) for (k,v) ∈ p.coeffs]...)


## -----
# unary subtraction
Base.:-(p::P) where {T,X,P<:FactoredPolynomial{T,X}} = (-1)*p

# addition
function Base.:+(p::P, c::S) where {S<:Number, T,X, P<:FactoredPolynomial{T,X}}
    R = promote_type(S,T)
    𝑷 = Polynomial{R,X}
    𝒑,𝒒 = convert(𝑷, p), convert(𝑷, c)
    convert(P, 𝒑 + 𝒒)
end
function Base.:+(p::P, q::P) where {T,X,P<:FactoredPolynomial{T,X}}
    𝑷 = Polynomial{T,X}
    𝒑,𝒒 = convert(𝑷, p), convert(𝑷, q)
    convert(P, 𝒑 + 𝒒 )
end

# subtraction
function Base.:-(p::P, c::S) where {S<:Number, T,X, P<:FactoredPolynomial{T,X}}
    R = promote_type(S,T)
    𝑷 = Polynomial{R,X}
    𝒑,𝒒 = convert(𝑷, p), convert(𝑷, c)
    convert(P, 𝒑 - 𝒒)
end
function Base.:-(p::P, q::P) where {T,X,P<:FactoredPolynomial{T,X}}
    𝑷 = Polynomial{T,X}
    𝒑,𝒒 = convert(𝑷, p), convert(𝑷, q)
    convert(P, 𝒑 - 𝒒 )
end

# multiplication
function Base.:*(p::P, q::P) where {T,X, P<:FactoredPolynomial{T,X}}
    d = copy(p.coeffs)
    for (k,v) ∈ q.coeffs
        d[k] = get(d, k, 0) + v
    end
    P(d, p.c*q.c)
end

# scalar mult
function scalar_mult(p::P, c::S) where {S<:Number, T, X, P <: FactoredPolynomial{T, X}}
    R = promote_type(T,S)
    d = OrderedDict{R, Int}() # wident
    for (k,v) ∈ p.coeffs
        d[k] = v
    end
    FactoredPolynomial{R,X}(d, c * p.c)
end
scalar_mult(c::S, p::P) where {S<:Number, T, X, P <: FactoredPolynomial{T, X}} = scalar_mult(p, c) # assume commutative, as we have S <: Number

# scalar division
function Base.:/(p::P, c::S) where {S<:Number, T, X, P <: FactoredPolynomial{T, X}}
    p * (1/c)
end

function Base.:^(p::P, n::Integer) where {T,X, P<:FactoredPolynomial{T,X}}
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    d = OrderedDict{T,Int}()
    for (k,v) ∈ p.coeffs
        d[k] = v*n
    end
    P(d, p.c^n)
end

## ----
## gcd, divrem, uvw
function Base.gcd(p::P, q::P) where {T, X, P<:FactoredPolynomial{T,X}}
    iszero(p) && return q
    iszero(q) && return p
    d = OrderedDict{T,Int}()

    for k ∈ intersect(keys(p.coeffs), keys(q.coeffs))
        d[k] = min(p.coeffs[k], q.coeffs[k])
    end

    P(d)
end

# return u,v,w with p = u*v , q = u*w
function uvw(p::P, q::P; kwargs...) where {T, X, P<:FactoredPolynomial{T,X}}
    du, dv, dw = OrderedDict{T,Int}(), OrderedDict{T,Int}(), OrderedDict{T,Int}()
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
    P(du), P(dv, p.c), P(dw, q.c)
end

# return a,b with p = a + b*q
function Base.divrem(p::P, q::P) where {T, X, P<:FactoredPolynomial{T,X}}

    u, v, w = uvw(p, q)
    𝑷 = Polynomial
    𝒗,𝒘 = convert(𝑷, v), convert(𝑷, w)
    𝒅,𝒓 = divrem(𝒗, 𝒘)
    d,r = convert(P, 𝒅), convert(P, 𝒓)

    return (d, u*r)

end

## ----

function integrate(p::P) where {P <: FactoredPolynomial}
    𝑷 = Polynomial
    𝒑 = convert(𝑷, p)
    ∫𝒑 = integrate(𝒑)
    convert(P, ∫𝒑)
end

function derivative(p::P,n::Int=1) where {P <: FactoredPolynomial}
    𝑷 = Polynomial
    𝒑 = convert(𝑷, p)
    𝒑⁽ⁿ⁾ = derivative(𝒑, n)
    convert(P, 𝒑⁽ⁿ⁾)
end

function Multroot.multroot(p::FactoredPolynomial)
    d = p.coeffs
    (values = collect(keys(d)),
     multiplicities = collect(values(d)),
     κ = 0.0, ϵ = 0.0)
end
