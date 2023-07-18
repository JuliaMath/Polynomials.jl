struct MutableSparsePolynomial{B,T,X} <:  AbstractUnivariatePolynomial{B, T,X}
    coeffs::Dict{Int, T}
    function MutableSparsePolynomial{B,T,X}(coeffs::AbstractDict{Int,T},order::Int=0) where {B,T,X}
        for (i, cᵢ) ∈ pairs(coeffs)
            iszero(cᵢ) && delete!(coeffs, i)
        end
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
function MutableSparsePolynomial{B,T,X}(coeffs::AbstractDict{Int,S}) where {B,T,S,X}
    cs = convert(Dict{Int,T}, coeffs)
    MutableSparsePolynomial{B,T,X}(cs)
end

function MutableSparsePolynomial{B}(cs::AbstractDict{Int,T}, var::SymbolLike=:x) where {B,T}
    MutableSparsePolynomial{B,T,Symbol(var)}(cs)
end

# abstract vector has order/symbol
function MutableSparsePolynomial{B,T,X}(coeffs::AbstractVector{S}, order::Int=0) where {B,T,S,X}

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

function MutableSparsePolynomial{B,T}(coeffs::AbstractVector{S}, var::SymbolLike) where {B,T,S}
    MutableSparsePolynomial{B,T,Symbol(var)}(xs)
end

function MutableSparsePolynomial{B}(xs::AbstractVector{T}, order::Int, var::SymbolLike=Var(:x)) where {B,T}
    MutableSparsePolynomial{B,T,Symbol(var)}(xs)
end

function MutableSparsePolynomial{B}(xs::AbstractVector{T}, var::SymbolLike) where {B,T}
    MutableSparsePolynomial{B,T,Symbol(var)}(xs)
end

function MutableSparsePolynomial{B,T}(xs, var::SymbolLike=Var(:x)) where {B,T}
    cs = collect(T, xs)
    cs = trim_trailing_zeros(cs)
    MutableSparsePolynomial{B,T,Symbol(var)}(cs)
end

function MutableSparsePolynomial{B}(xs, var::SymbolLike=Var(:x)) where {B}
    cs = collect(xs)
    cs = trim_trailing_zeros(cs)
    MutableSparsePolynomial{B,T,Symbol(var)}(cs)
end


@poly_register MutableSparsePolynomial
constructorof(::Type{<:MutableSparsePolynomial}) = MutableSparsePolynomial


# cs iterable of pairs
# XXX enure tight value of T
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
        for (a,b) ∈ c
            T = promote_type(T, typeof(b))
        end
        ks = Base.Generator(first, cs)
        vs = Base.Generator(last, cs)
        d = Dict{Int,T}(Base.Generator(=>, ks, vs))
        return MutableSparsePolynomial{B,T,X}(d)
    end
end



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
Base.iterate(p::MutableSparsePolynomial, args...) = Base.iterate(p.coeffs, args...)
Base.pairs(p::MutableSparsePolynomial) = pairs(p.coeffs)

## Not properly named!!! truncate? chop? skip?
function trim_trailing_zeros(d::Dict)
    for (k,v) ∈ pairs(d)
        iszero(v) && deletat!(d, k)
    end
    d
end

function Base.chop(p::MutableSparsePolynomial; atol=nothing, rtol=nothing)
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, _norm(x,2) * δ)
    for (i,pᵢ) ∈ pairs(p)
        abs(pᵢ) ≤ τ && delete!(p.coeffs, i)
    end
    p
end

_zeros(::Type{MutableSparsePolynomial{B,T,X}}, z::S, N) where {B,T,X,S} = Dict{Int, S}()

Base.iszero(p::MutableSparsePolynomial) = all(iszero, values(p.coeffs))

Base.zero(::Type{MutableSparsePolynomial{B,T,X}}) where {B,T,X} = MutableSparsePolynomial{B,T,X}(Dict{Int,T}())
## ---

function _evalpoly(p::MutableSparsePolynomial, x)

    tot = zero(p[0]*x)
    for (i, cᵢ) ∈ p.coeffs
        tot = muladd(cᵢ, x^i, tot)
    end
    return tot
end

offset(p::MutableSparsePolynomial) = 0
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

function scalar_mul(c::S, p::MutableSparsePolynomial{B,T,X}) where {B,T,X,S}

    R = promote_type(T,S)
    P = MutableSparsePolynomial{B,R,X}
    (iszero(p) || iszero(c)) && return(zero(P))

    d = convert(Dict{Int, R}, copy(p.coeffs))
    for (k, pₖ) ∈ pairs(d)
        @inbounds d[k] = c .* d[k]
    end

    return P(Val(false), d)

end

function scalar_mul(p::MutableSparsePolynomial{B,T,X}, c::S) where {B,T,X,S}
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
