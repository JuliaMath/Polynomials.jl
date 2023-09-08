"""
    MutableSparseVectorPolynomial{B,T,X}

This polynomial type uses an `SparseVector{T,Int}` to store the coefficients of a polynomial relative to the basis `B` with indeterminate `X`.
The type `T` should have `zero(T)` defined.


"""
struct MutableSparseVectorPolynomial{B,T,X} <:  AbstractUnivariatePolynomial{B, T,X}
    coeffs::SparseVector{T, Int}
    function MutableSparseVectorPolynomial{B,T,X}(cs::SparseVector{S,Int}, order::Int=0) where {B,T,S,X}
        new{B,T,Symbol(X)}(cs)
    end
end

MutableSparseVectorPolynomial{B,T,X}(check::Val{:false}, coeffs::SparseVector{Int,S}) where {B,T,S,X} =
    MutableSparseVectorPolynomial{B,T,X}(coeffs)
MutableSparseVectorPolynomial{B,T,X}(checked::Val{:true}, coeffs::SparseVector{Int,T}) where {B,T,X<:Symbol} =
    MutableSparseVectorPolynomial{B,T,X}(coeffs)

# ---
function MutableSparseVectorPolynomial{B,T}(coeffs::SparseVector{S,Int}, var::SymbolLike=Var(:x)) where {B,T,S}
    MutableSparseVectorPolynomial{B,T,Symbol(var)}(coeffs)
end

function MutableSparseVectorPolynomial{B}(cs::SparseVector{T,Int}, var::SymbolLike=Var(:x)) where {B,T}
    MutableSparseVectorPolynomial{B,T,Symbol(var)}(cs)
end

# From a Dictionary
function MutableSparseVectorPolynomial{B,X}(cs::AbstractDict{Int, T}) where {B,T,X}
    N = maximum(keys(cs)) + 1
    v = SparseVector(N, 1 .+ keys(cs), collect(values(cs)))
    MutableSparseVectorPolynomial{B,T,X}(v)
end

function MutableSparseVectorPolynomial{B}(cs::AbstractDict{Int, T}, var::SymbolLike=Var(:x)) where {B,T}
    MutableSparseVectorPolynomial{B,Symbol(var)}(cs)
end


# abstract vector has order/symbol
function MutableSparseVectorPolynomial{B,T,X}(coeffs::AbstractVector{S}, order::Int=0) where {B,T,S,X}
    if Base.has_offset_axes(coeffs)
        @warn "ignoring the axis offset of the coefficient vector"
        coeffs = parent(coeffs)
    end

    MutableSparseVectorPolynomial{B,T,X}(convert(SparseVector, coeffs))
end


# # cs iterable of pairs; ensuring tight value of T
# function MutableSparseVectorPolynomial{B}(cs::Tuple, var::SymbolLike=:x) where {B}
#     isempty(cs) && throw(ArgumentError("No type attached"))
#     X = Var(var)
#     if length(cs) == 1
#         c = only(cs)
#         d = Dict(first(c) => last(c))
#         T = eltype(last(c))
#         return MutableSparseVectorPolynomial{B,T,X}(d)
#     else
#         c₁, c... = cs
#         T = typeof(last(c₁))
#         for b ∈ c
#             T = promote_type(T, typeof(b))
#         end
#         ks = 0:length(cs)-1
#         vs = cs
#         d = Dict{Int,T}(Base.Generator(=>, ks, vs))
#         return MutableSparseVectorPolynomial{B,T,X}(d)
#     end
# end

constructorof(::Type{<:MutableSparseVectorPolynomial{B}}) where {B <: AbstractBasis} = MutableSparseVectorPolynomial{B}
@poly_register MutableSparseVectorPolynomial

function Base.map(fn, p::P, args...)  where {B,T,X, P<:MutableSparseVectorPolynomial{B,T,X}}
    xs = map(fn, p.coeffs)
    R = eltype(xs)
    return MutableSparseVectorPolynomial{B, R, X}(xs)
end

function Base.map!(fn, q::Q, p::P, args...)  where {B,T,X, P<:MutableSparseVectorPolynomial{B,T,X},S,Q<:MutableSparseVectorPolynomial{B,S,X}}
    map!(fn, p.coeffs, p.coeffs)
    nothing
end

## ---
Base.collect(p::MutableSparseVectorPolynomial) = collect(p.coeffs)
Base.collect(::Type{T}, p::MutableSparseVectorPolynomial) where {T} = collect(T, p.coeffs)
minimumexponent(::Type{<:MutableSparseVectorPolynomial}) =  0

Base.length(p::MutableSparseVectorPolynomial) = length(p.coeffs)

function degree(p::MutableSparseVectorPolynomial)
    idx = findall(!iszero, p.coeffs)
    isempty(idx) && return -1
    n = maximum(idx)
    n - 1
end

Base.copy(p::MutableSparseVectorPolynomial{B,T,X}) where {B,T,X} = MutableSparseVectorPolynomial{B,T,X}(copy(p.coeffs))

function Base.convert(::Type{MutableSparseVectorPolynomial{B,T,X}}, p::MutableSparseVectorPolynomial{B,S,X}) where {B,T,S,X}
    cs = convert(SparseVector{T,Int}, p.coeffs)
    MutableSparseVectorPolynomial{B,T,X}(cs)
end

function Base.:(==)(p1::P, p2::P) where {P <: MutableSparseVectorPolynomial}
    iszero(p1) && iszero(p2) && return true

    ks1 = findall(!iszero, p1.coeffs)
    ks2 = findall(!iszero, p2.coeffs)
    length(ks1) == length(ks2) || return false
    idx = sortperm(ks1)
    for i ∈ idx
        ks1[i] == ks2[i] || return false
        p1.coeffs[ks1[i]] == p2.coeffs[ks2[i]] || return false
    end

    return true
    # #    eachindex(p1) == eachindex(p2) || return false
    # # coeffs(p1) == coeffs(p2), but non-allocating
    # p1val = (p1[i] for i in eachindex(p1))
    # p2val = (p2[i] for i in eachindex(p2))
    # all(((a,b),) -> a == b, zip(p1val, p2val))
end

# ---

Base.firstindex(p::MutableSparseVectorPolynomial) = 0
function Base.lastindex(p::MutableSparseVectorPolynomial)
    isempty(p.coeffs) && return 0
    maximum(keys(p.coeffs))
end

function Base.getindex(p::MutableSparseVectorPolynomial{B,T,X}, i::Int) where {B,T,X}
    get(p.coeffs, i + 1, zero(T))
end

# errors if extending
function Base.setindex!(p::MutableSparseVectorPolynomial{B,T,X}, value, i::Int) where {B,T,X}
    p.coeffs[i+1] = value
end


function Base.pairs(p::MutableSparseVectorPolynomial)
    ks, vs = findnz(p.coeffs)
    idx = sortperm(ks) # guarantee order here
    Base.Generator(=>, ks[idx] .- 1, vs)
end
Base.keys(p::MutableSparseVectorPolynomial) = Base.Generator(first, pairs(p))
Base.values(p::MutableSparseVectorPolynomial) = Base.Generator(last, pairs(p))

basis(P::Type{<:MutableSparseVectorPolynomial{B, T, X}}, i::Int) where {B,T,X} = P(SparseVector(1+i, [i+1], [1]))

# return coeffs as  a vector
function coeffs(p::MutableSparseVectorPolynomial{B,T})  where {B,T}
    d = degree(p)
    ps = p.coeffs
    [ps[i] for i ∈ 1:(d+1)]
end


hasnan(p::MutableSparseVectorPolynomial) = any(hasnan, values(p.coeffs))::Bool


offset(p::MutableSparseVectorPolynomial) = 1

function keys_union(p::MutableSparseVectorPolynomial, q::MutableSparseVectorPolynomial)
    # IterTools.distinct(Base.Iterators.flatten((keys(p), keys(q)))) may allocate less
    unique(Base.Iterators.flatten((keys(p), keys(q))))
end



## ---

chop_exact_zeros!(d::SparseVector{T, Int}) where {T} = d


function _truncate!(v::SparseVector{T,X};
                    rtol::Real = Base.rtoldefault(real(T)),
                    atol::Real = 0) where {T,X}
    isempty(v) && return v
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, norm(values(v),2) * δ)
    for (i,pᵢ) ∈ pairs(v)
        abs(pᵢ) ≤ τ && (v[i] = zero(T))
    end
    v
end


chop!(p::MutableSparseVectorPolynomial; kwargs...) = (chop!(p.coeffs; kwargs...); p)
function chop!(d::SparseVector{T, Int}; atol=nothing, rtol=nothing) where {T}
    isempty(d) && return d
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, norm(values(d),2) * δ)
    for (i, pᵢ) ∈ Base.Iterators.reverse(pairs(d))
        abs(pᵢ) ≥ τ && break
        d[i] = zero(T)
    end
    d
end

## ---

_zeros(::Type{MutableSparseVectorPolynomial{B,T,X}}, z::S, N) where {B,T,X,S} = zeros(T, N)

Base.zero(::Type{MutableSparseVectorPolynomial{B,T,X}}) where {B,T,X} = MutableSparseVectorPolynomial{B,T,X}(spzeros(T,0))

## ---

function isconstant(p::MutableSparseVectorPolynomial)
    degree(p) <= 0
end

Base.:+(p::MutableSparseVectorPolynomial{B,T,X}, q::MutableSparseVectorPolynomial{B,S,X}) where{B,X,T,S} =
    _sparse_vector_combine(+, p, q)
Base.:-(p::MutableSparseVectorPolynomial{B,T,X}, q::MutableSparseVectorPolynomial{B,S,X}) where{B,X,T,S} =
    _sparse_vector_combine(-, p, q)

# embed into bigger vector
function _embed(v::SparseVector{T, Int}, l) where {T}
    l == length(v) && return v
    ks,vs = findnz(v)
    SparseVector(l, ks, vs)
end


function _sparse_vector_combine(op, p::MutableSparseVectorPolynomial{B,T,X}, q::MutableSparseVectorPolynomial{B,S,X}) where{B,X,T,S}
    R = promote_type(T,S)
    ps, qs = p.coeffs, q.coeffs
    m = max(length(ps), length(qs))
    ps′, qs′ = _embed(ps, m), _embed(qs, m)
    cs = op(ps′, qs′)
    MutableSparseVectorPolynomial{B,R,X}(cs)
end
