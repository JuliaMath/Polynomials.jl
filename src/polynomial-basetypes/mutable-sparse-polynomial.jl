"""

This polynomial type uses an `Dict{Int,T}` to store the coefficients of a polynomial relative to the basis `B` with indeterminate `X`.
Explicit `0` coefficients are not stored. This type can be used for Laurent polynomials.

"""
struct MutableSparsePolynomial{B,T,X} <:  AbstractLaurentUnivariatePolynomial{B, T,X}
     coeffs::Dict{Int, T}
    function MutableSparsePolynomial{B,T,X}(cs::AbstractDict{Int,S},order::Int=0) where {B,T,S,X}
        coeffs = convert(Dict{Int,T}, cs)
        chop_exact_zeros!(coeffs)
        new{B,T,Symbol(X)}(coeffs)
    end
    function MutableSparsePolynomial{B,T,X}(check::Val{:false}, coeffs::AbstractDict{Int,S}) where {B,T,S,X}
        new{B,T,Symbol(X)}(coeffs)
    end
end

function MutableSparsePolynomial{B,T,X}(checked::Val{:true}, coeffs::AbstractDict{Int,T}) where {B,T,X<:Symbol}
    MutableSparsePolynomial{B,T,X}(coeffs)
end

function MutableSparsePolynomial{B,T}(coeffs::AbstractDict{Int,S}, var::SymbolLike=Var(:x)) where {B,T,S}
    MutableSparsePolynomial{B,T,Symbol(var)}(coeffs)
end

function MutableSparsePolynomial{B}(cs::AbstractDict{Int,T}, var::SymbolLike=Var(:x)) where {B,T}
    MutableSparsePolynomial{B,T,Symbol(var)}(cs)
end

# abstract vector has order/symbol
function MutableSparsePolynomial{B,T,X}(coeffs::AbstractVector{S}, order::Int=0) where {B,T,S,X}
    if Base.has_offset_axes(coeffs)
        @warn "ignoring the axis offset of the coefficient vector"
        coeffs = parent(coeffs)
    end

    P = MutableSparsePolynomial{B,T,X}
    n = length(coeffs)
    iszero(n) && zero(P)
    xs = convert(Vector{T}, coeffs)
    d = Dict{Int, T}(Base.Generator(=>, order:(order+n-1), xs))
    P(d)
end


# cs iterable of pairs; ensuring tight value of T
function MutableSparsePolynomial{B}(cs::Tuple, var::SymbolLike=:x) where {B}
    isempty(cs) && throw(ArgumentError("No type attached"))
    X = Var(var)
    if length(cs) == 1
        c = only(cs)
        d = Dict(first(c) => last(c))
        T = eltype(last(c))
        return MutableSparsePolynomial{B,T,X}(d)
    else
        c₁, c... = cs
        T = typeof(last(c₁))
        for b ∈ c
            T = promote_type(T, typeof(b))
        end
        ks = 0:length(cs)-1
        vs = cs
        d = Dict{Int,T}(Base.Generator(=>, ks, vs))
        return MutableSparsePolynomial{B,T,X}(d)
    end
end

@poly_register MutableSparsePolynomial

constructorof(::Type{<:MutableSparsePolynomial{B}}) where {B} = MutableSparsePolynomial{B}

function Base.map(fn, p::P, args...)  where {B,T,X, P<:MutableSparsePolynomial{B,T,X}}
    xs = Dict(k => fn(v, args...) for (k,v) ∈ pairs(p.coeffs))
    xs = chop_exact_zeros!(xs)
    R = eltype(values(xs)) # narrow_eltype...
    return ⟒(P){R,X}(Val(false), xs)
end

## ---

minimumexponent(::Type{<:MutableSparsePolynomial}) =  typemin(Int)

Base.copy(p::MutableSparsePolynomial{B,T,X}) where {B,T,X} = MutableSparsePolynomial{B,T,X}(copy(p.coeffs))

function Base.convert(::Type{MutableSparsePolynomial{B,T,X}}, p::MutableSparsePolynomial{B,S,X}) where {B,T,S,X}
    d = Dict{Int,T}(k => v for (k,v) ∈ pairs(p.coeffs))
    MutableSparsePolynomial{B,T,X}(Val(false), d)
end

# ---

function Base.firstindex(p::MutableSparsePolynomial)
    isempty(p.coeffs) && return 0
    i = minimum(keys(p.coeffs))
end

function Base.lastindex(p::MutableSparsePolynomial)
    isempty(p.coeffs) && return 0
    maximum(keys(p.coeffs))
end

function Base.getindex(p::MutableSparsePolynomial{B,T,X}, i::Int) where {B,T,X}
    get(p.coeffs, i, zero(T))
end

function Base.setindex!(p::MutableSparsePolynomial{B,T,X}, value, i::Int) where {B,T,X}
    iszero(value) && delete!(p.coeffs, i)
    p.coeffs[i] = value
end

# would like this, but fails a test... (iterate does not guarantee any order)
#Base.iterate(p::MutableSparsePolynomial, args...) = throw(ArgumentError("Use `pairs` to iterate a sparse polynomial"))

# return coeffs as  a vector
# use p.coeffs to get Dictionary
function coeffs(p::MutableSparsePolynomial{B,T})  where {B,T}
    a,b = min(0,firstindex(p)), lastindex(p)
    cs = zeros(T, length(a:b))
    for k in sort(collect(keys(p.coeffs)))
        v = p.coeffs[k]
        cs[k - a + 1] = v
    end
    cs
end


hasnan(p::MutableSparsePolynomial) = any(hasnan, values(p.coeffs))
Base.pairs(p::MutableSparsePolynomial) = pairs(p.coeffs)

offset(p::MutableSparsePolynomial) = 0

## ---

function chop_exact_zeros!(d::Dict)
    for (k,v) ∈ pairs(d)
        iszero(v) && delete!(d, k)
    end
    d
end
trim_trailing_zeros!!(d::Dict) = chop_exact_zeros!(d) # Not properly named, but what is expected in other constructors

chop!(p::MutableSparsePolynomial; kwargs...) = (chop!(p.coeffs; kwargs...); p)
function chop!(p::Dict; atol=nothing, rtol=nothing)
    isempty(p.coeffs) && return p
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, norm(values(p.coeffs),2) * δ)
    for (i,pᵢ) ∈ pairs(p)
        abs(pᵢ) ≤ τ && delete!(p.coeffs, i)
    end
    p
end

## ---

_zeros(::Type{MutableSparsePolynomial{B,T,X}}, z::S, N) where {B,T,X,S} = Dict{Int, S}()

Base.zero(::Type{MutableSparsePolynomial{B,T,X}}) where {B,T,X} = MutableSparsePolynomial{B,T,X}(Dict{Int,T}())

## ---

function isconstant(p::MutableSparsePolynomial)
    n = length(p.coeffs)
    n == 0 && return true
    n == 1 && haskey(p.coeffs, 0)
end

Base.:+(p::MutableSparsePolynomial{B,T,X}, q::MutableSparsePolynomial{B,S,X}) where{B,X,T,S} =
    _dict_combine(+, p, q)
Base.:-(p::MutableSparsePolynomial{B,T,X}, q::MutableSparsePolynomial{B,S,X}) where{B,X,T,S} =
    _dict_combine(-, p, q)

function _dict_combine(op, p::MutableSparsePolynomial{B,T,X}, q::MutableSparsePolynomial{B,S,X}) where{B,X,T,S}

    R = promote_type(T,S)
    P = MutableSparsePolynomial{B,R,X}
    D = convert(Dict{Int, R}, copy(p.coeffs))
    for (i, qᵢ) ∈ pairs(q.coeffs)
        pᵢ =  get(D, i, zero(R))
        pqᵢ = op(pᵢ, qᵢ)
        if iszero(pqᵢ)
            delete!(D, i) # will be zero
        else
            D[i] = pqᵢ
        end
    end
    return P(Val(false), D)

end
