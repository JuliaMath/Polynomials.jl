# Try to keep length based on N,M so no removal of trailing zeros by default
# order is ignored, firstindex is always 0
struct ImmutableDensePolynomial{B,T,X,N} <: AbstractUnivariatePolynomial{B,T,X}
    coeffs::NTuple{N,T}
    function ImmutableDensePolynomial{B,T,X,N}(cs::NTuple{N,T}) where {B,N,T,X}
        if Base.has_offset_axes(cs)
            @warn "Ignoring the axis offset of the coefficient vector"
        end
        new{B,T,Symbol(X),N}(cs)
    end
end

ImmutableDensePolynomial{B,T,X,N}(::Type{Val{false}}, cs::NTuple{N,T}) where {B,N,T,X} =
    ImmutableDensePolynomial{B,T,X}(cs)

ImmutableDensePolynomial{B,T,X,N}(::Type{Val{true}}, cs::NTuple{N,T}) where {B,N, T,X} =
    ImmutableDensePolynomial{B,T,X,N}(cs)

# tuple with mis-matched size
function ImmutableDensePolynomial{B,T,X,N}(xs::NTuple{M,S}) where {B,T,S,X,N,M}
    convert(ImmutableDensePolynomial{B,T,X,N}, ImmutableDensePolynomial{B,T,X,M}(xs))
end

function ImmutableDensePolynomial{B,T,X}(xs::NTuple{N,S}) where {B,T,S,X,N}
    cs = convert(NTuple{N,T}, xs)
    ImmutableDensePolynomial{B,T,X,N}(cs)
end

function ImmutableDensePolynomial{B,T}(xs::NTuple{N,S}, var::SymbolLike=Var(:x)) where {B,T,S,N}
    ImmutableDensePolynomial{B,T,Var(var),N}(xs)
end

function ImmutableDensePolynomial{B}(xs::NTuple{N,T}, var::SymbolLike=Var(:x)) where {B,T,N}
    ImmutableDensePolynomial{B,T,Var(var),N}(xs)
end

# abstract vector
function ImmutableDensePolynomial{B,T,X}(xs::AbstractVector{S}, order::Int=0) where {B,T,X,S}
    if Base.has_offset_axes(xs)
        @warn "ignoring the axis offset of the coefficient vector"
        xs = parent(xs)
    end
    !iszero(order) && @warn "order argument is ignored"
    N = length(xs)
    cs = ntuple(Base.Fix1(T∘getindex,xs), Val(N))
    ImmutableDensePolynomial{B,T,X,N}(cs)
end

function ImmutableDensePolynomial{B,T}(xs::AbstractVector{S}, order::Int, var::SymbolLike=Var(:x)) where {B,T,S}
    ImmutableDensePolynomial{B,T,Symbol(var),N}(xs)
end

function ImmutableDensePolynomial{B,T}(xs::AbstractVector{S}, var::SymbolLike) where {B,T,S}
    ImmutableDensePolynomial{B,T,Symbol(var)}(xs)
end

function ImmutableDensePolynomial{B}(xs::AbstractVector{T}, order::Int, var::SymbolLike=Var(:x)) where {B,T}
    ImmutableDensePolynomial{B,T,Symbol(var)}(xs)
end

function ImmutableDensePolynomial{B}(xs::AbstractVector{T}, var::SymbolLike=Var(:x)) where {B,T}
    ImmutableDensePolynomial{B,T,Symbol(var)}(xs)
end

function ImmutableDensePolynomial{B,T,X}(xs) where {B,T,X}
    cs = collect(T,xs)
    ImmutableDensePolynomial{B,T,X}(cs)
end

function ImmutableDensePolynomial{B,T}(xs; var::SymbolLike=Var(:x)) where {B,T}
    cs = collect(T,xs)
    ImmutableDensePolynomial{B,T,Var(var)}(cs)
end


function ImmutableDensePolynomial{B}(xs, var::SymbolLike=Var(:x)) where {B}
    cs = collect(promote(xs...))
    T = eltype(cs)
    ImmutableDensePolynomial{B,T,Symbol(var)}(cs)
end

# constant
function ImmutableDensePolynomial{B,T,X,N}(c::S) where {B,T,X,N,S<:Number}
    if N == 0
        if iszero(c)
            throw(ArgumentError("Can't create zero-length polynomial"))
        else
            return zero(ImmutableDensePolynomial{B,T,X})
        end
    end
    cs = ntuple(i -> i == 1 ? T(c) : zero(T), Val(N))
    return ImmutableDensePolynomial{B,T,X,N}(cs)
end

@poly_register ImmutableDensePolynomial
constructorof(::Type{<:ImmutableDensePolynomial{B}})  where {B} = ImmutableDensePolynomial{B}

## ----

# need to promote to larger
Base.promote_rule(::Type{<:ImmutableDensePolynomial{B,T,X,N}}, ::Type{<:ImmutableDensePolynomial{B,S,X,M}}) where {B,T,S,X,N,M} =
    ImmutableDensePolynomial{B,promote_type(T,S), X, max(N,M)}
Base.promote_rule(::Type{<:ImmutableDensePolynomial{B,T,X,N}}, ::Type{<:S}) where {B,T,S<:Number,X,N} =
    ImmutableDensePolynomial{B,promote_type(T,S), X, N}


function Base.convert(::Type{<:ImmutableDensePolynomial{B,T,X,N}},
                      p::ImmutableDensePolynomial{B,T′,X,N′}) where {B,T,T′,X,N,N′}
    N′′ = findlast(!iszero, p)
    isnothing(N′′) && return zero(ImmutableDensePolynomial{B,T,X,N})
    N < N′′ && throw(ArgumentError("Wrong size"))
    N > N′′ && return ImmutableDensePolynomial{B,T,X,N}(ntuple(i -> T(p[i-1]), Val(N)))
    ImmutableDensePolynomial{B,T,X,N}(ntuple(i -> T(p[i-1]), Val(N)))
end


Base.copy(p::ImmutableDensePolynomial) = p
Base.similar(p::ImmutableDensePolynomial, args...) = p.coeffs

## chop
function Base.chop(p::ImmutableDensePolynomial{B,T,X,N};
                   rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0) where {B,T,X,N}
    i = chop_right_index(p.coeffs; rtol=rtol, atol=atol)
    if i == nothing
        xs = ()
        N′ = 0
    else
        N′ = i
        xs = ntuple(Base.Fix1(getindex, p.coeffs), Val(N′))
    end
    ImmutableDensePolynomial{B,T,X,N′}(xs)
end

# misnamed!
chop!(p::ImmutableDensePolynomial; kwargs...) = chop(p; kwargs...)

function truncate!(p::ImmutableDensePolynomial{B,T,X,N};
                   rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0) where {B,T,X,N}

    ps = p.coeffs
    isempty(ps) && return p
    max_coeff = norm(ps, Inf)
    thresh = max_coeff * rtol + atol
    for (i,pᵢ) ∈ pairs(ps)
        if abs(pᵢ) <= thresh
            @set! ps[i] = zero(T)
        end
    end
    ImmutableDensePolynomial{B,T,X,N}(ps)
end


# not type stable, as N is value dependent
function trim_trailing_zeros(cs::Tuple)
    isempty(cs) && return cs
    !iszero(last(cs)) && return cs
    i = findlast(!iszero, cs)
    i == nothing && return ()
    xs = ntuple(Base.Fix1(getindex,cs), Val(i))
    xs
end

# isapprox helper
function normΔ(q1::ImmutableDensePolynomial{B}, q2::ImmutableDensePolynomial{B}, p::Real = 2) where {B}
    iszero(q1) && return norm(q2, p)
    iszero(q2) && return norm(q1, p)
    r = zero(q1[end] + q2[end])
    tot = zero(r)
    for i ∈ 1:maximum(lastindex, (q1,q2))
       @inbounds tot += (q1[i] - q2[i])^p
    end
    return tot^(1/p)
end

# function Base.isapprox(p1::ImmutableDensePolynomial{B,T,X}, p2::ImmutableDensePolynomial{B,T′,X};
#                        atol=nothing, rtol = nothing
#                        ) where {B,T,T′,X}

#     (hasnan(p1) || hasnan(p2)) && return false  # NaN poisons comparisons
#     R = real(float(promote_type(T,T′)))
#     atol = something(atol, zero(R))
#     rtol = something(rtol, Base.rtoldefault(R))
#     # copy over from abstractarray.jl
#     Δ  = normΔ(p1,p2)
#     if isfinite(Δ)
#         return Δ <= max(atol, rtol * max(norm(p1), norm(p2)))
#     else
#         for i in 0:max(degree(p1), degree(p2))
#             isapprox(p1[i], p2[i]; rtol=rtol, atol=atol) || return false
#         end
#         return true
#     end
# end

## ---

_zeros(::Type{<:ImmutableDensePolynomial}, z::S, N) where {S} =
    ntuple(_ -> zero(S), Val(N))


Base.firstindex(p::ImmutableDensePolynomial) = 0
Base.lastindex(p::ImmutableDensePolynomial{B,T,X,N}) where {B,T,X,N} = N - 1
Base.eachindex(p::ImmutableDensePolynomial) = firstindex(p):lastindex(p)
Base.pairs(p::ImmutableDensePolynomial) =
    Base.Generator(=>, eachindex(p), p.coeffs)
Base.length(p::ImmutableDensePolynomial{B,T,X,N}) where {B,T,X,N} = N
offset(p::ImmutableDensePolynomial) = 1

function Base.getindex(p::ImmutableDensePolynomial{B,T,X,N}, i::Int) where {B,T,X,N}
    N == 0 && return zero(T)
    (i < firstindex(p) || i > lastindex(p)) && return zero(p.coeffs[1])
    p.coeffs[i + offset(p)]
end

Base.setindex!(p::ImmutableDensePolynomial, value, i::Int) =
    throw(ArgumentError("ImmutableDensePolynomial has no setindex! method"))

# can't promote to same N if trailing zeros
function Base.:(==)(p1::ImmutableDensePolynomial{B}, p2::ImmutableDensePolynomial{B}) where {B}
    iszero(p1) && iszero(p2) && return true
    if isconstant(p1)
        isconstant(p2) && return constantterm(p1) == constantterm(p2)
        return false
    elseif isconstant(p2)
        return false # p1 is not constant
    end
    check_same_variable(p1, p2) || return false
    for i ∈ union(eachindex(p1), eachindex(p2))
        p1[i] == p2[i] || return false
    end
    return true
end

minimumexponent(::Type{<:ImmutableDensePolynomial}) =  0


## ---
function degree(p::ImmutableDensePolynomial{B,T,X,N}) where {B,T,X,N}
    iszero(N) && return -1
    i = findlast(!iszero, p.coeffs)
    isnothing(i) && return -1
    return i - 1
end

# zero, one
Base.zero(::Type{<:ImmutableDensePolynomial{B,T,X}}) where {B,T,X} =
    ImmutableDensePolynomial{B,T,X,0}(())

Base.zero(::Type{ImmutableDensePolynomial{B,T,X,N}}) where {B,T,X,N} =
    ImmutableDensePolynomial{B,T,X,0}(())

function basis(P::Type{<:ImmutableDensePolynomial{B,T,X}}, i::Int) where {B,T,X}
    xs = zeros(T, i + 1)
    xs[end] = 1
    ImmutableDensePolynomial{B,T,X}(xs)
end

coeffs(p::ImmutableDensePolynomial) = p.coeffs



## Vector space operations
# vector ops +, -, c*x
## unary
Base.:-(p::ImmutableDensePolynomial{B,T,X,N}) where {B,T,X,N} =
    ImmutableDensePolynomial{B,T,X,N}(map(-, p.coeffs))

## binary
function Base.:+(p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where{B,X,T,S,N,M}
    N < M && return q + p
    _tuple_combine(+, p, q)
end
function Base.:-(p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where{B,X,T,S,N,M}
    N < M && return (-q) + p
    _tuple_combine(-, p, q)
end

# handle +, -; Assum N >= M
function _tuple_combine(op, p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where{B,X,T,S,N,M}

    @assert N >= M
    R = promote_type(T,S)
    P = ImmutableDensePolynomial{B,R,X}

    iszero(p) && return zero(P)
    #xs = ntuple(i -> i <= M ? R(op(p.coeffs[i],q.coeffs[i])) : R(p.coeffs[i]), Val(N))
    xs = _tuple_combine(op, p.coeffs, q.coeffs)
    P{N}(xs)

end

# Padded vector combination of two homogeneous tuples assuming N ≥ M
@generated function _tuple_combine(op, p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}

    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(op(p1[$i],p2[$i]))
    end
    for i in M+1:N
        exprs[i] =:(p1[$i])
    end

    return quote
        Base.@_inline_meta
        #Base.@inline
        tuple($(exprs...))
    end

end

# scalar

function scalar_mult(p::ImmutableDensePolynomial{B,T,X,N}, c::S) where {B,T,X,S,N}
    iszero(N) && return zero(ImmutableDensePolynomial{B,T,X})
    iszero(c) && ImmutableDensePolynomial{B}([p[0] .* c], X)
    cs = p.coeffs .* (c,)
    R = eltype(cs)
    return ImmutableDensePolynomial{B,R,X,N}(cs)
end

function scalar_mult(c::S, p::ImmutableDensePolynomial{B,T,X,N}) where {B,T,X,S,N}
    iszero(N) && return zero(ImmutableDensePolynomial{B,T,X})
    iszero(c) && ImmutableDensePolynomial{B}([c .* p[0]],X)
    cs = (c,) .* p.coeffs
    R = eltype(cs)
    return ImmutableDensePolynomial{B,R,X,N}(cs)
end
