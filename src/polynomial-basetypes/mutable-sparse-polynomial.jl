# dictionary to store (i, cᵢ)
# ensure cᵢ ≠ 0 in constructor
struct MutableSparsePolynomial{B,T,X} <:  AbstractUnivariatePolynomial{B, T,X}
    coeffs::Dict{Int, T}
    function MutableSparsePolynomial{B,T,X}(cs::AbstractDict{Int,S},order::Int=0) where {B,T,S,X}
        coeffs = convert(Dict{Int,T}, cs)
        chop_exact_zeros!(coeffs)
        new{B,T,Symbol(X)}(coeffs)
    end
    function MutableSparsePolynomial{B,T,X}(checked::Val{:false}, coeffs::AbstractDict{Int,T},order::Int=0) where {B,T,X}
        new{B,T,Symbol(X)}(coeffs)
    end
end

function MutableSparsePolynomial{B,T,X}(checked::Val{:true}, coeffs::AbstractDict{Int,T}) where {B,T,X<:Symbol}
    MutableSparsePolynomial{B,T,X}(coeffs)
end

# Dict
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

function MutableSparsePolynomial{B,T}(xs::AbstractVector{S}, order::Int, var::SymbolLike=Var(:x)) where {B,T,S}
    MutableSparsePolynomial{B,T,Symbol(var)}(xs)
end

function MutableSparsePolynomial{B,T}(xs::AbstractVector{S}, var::SymbolLike) where {B,T,S}
    MutableSparsePolynomial{B,T,Symbol(var)}(xs)
end

function MutableSparsePolynomial{B}(xs::AbstractVector{T}, order::Int, var::SymbolLike=Var(:x)) where {B,T}
    MutableSparsePolynomial{B,T,Symbol(var)}(xs)
end

function MutableSparsePolynomial{B}(xs::AbstractVector{T}, var::SymbolLike) where {B,T}
    MutableSparsePolynomial{B,T,Symbol(var)}(xs)
end

# iterable
function MutableSparsePolynomial{B,T}(xs, var::SymbolLike=Var(:x)) where {B,T}
    cs = collect(T, xs)
    cs = trim_trailing_zeros(cs)
    MutableSparsePolynomial{B,T,Symbol(var)}(cs)
end

function MutableSparsePolynomial{B}(xs, var::SymbolLike=Var(:x)) where {B}
    cs = collect(xs)
    cs = trim_trailing_zeros(cs)
    MutableSparsePolynomial{B,eltype(cs),Symbol(var)}(cs)
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
trim_trailing_zeros(d::Dict) = chop_exact_zeros!(d) # Not properly named, but what is expected in other constructors

function chop!(p::MutableSparsePolynomial; atol=nothing, rtol=nothing)
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

# much faster than default
function scalar_add(c::S, p::MutableSparsePolynomial{B,T,X}) where {B,T,X,S}
    c₀ = c + p[0]
    R = eltype(c₀)
    P = MutableSparsePolynomial{B,R,X}
    D = convert(Dict{Int, R}, copy(p.coeffs))
    if iszero(c₀)
        delete!(D,0)
    else
        @inbounds D[0] = c₀
    end
    return P(Val(false), D)
end

function scalar_mult(c::S, p::MutableSparsePolynomial{B,T,X}) where {B,T,X,S}

    R = promote_type(T,S)
    P = MutableSparsePolynomial{B,R,X}
    (iszero(p) || iszero(c)) && return(zero(P))

    d = convert(Dict{Int, R}, copy(p.coeffs))
    for (k, pₖ) ∈ pairs(d)
        @inbounds d[k] = c .* d[k]
    end

    return P(Val(false), d)

end

function scalar_mult(p::MutableSparsePolynomial{B,T,X}, c::S) where {B,T,X,S}
        R = promote_type(T,S)
    P = MutableSparsePolynomial{B,R,X}
    (iszero(p) || iszero(c)) && return(zero(P))

    d = convert(Dict{Int, R}, copy(p.coeffs))
    for (k, pₖ) ∈ pairs(d)
        @inbounds d[k] = d[k] .* c
    end

    return P(Val(false), d)
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
