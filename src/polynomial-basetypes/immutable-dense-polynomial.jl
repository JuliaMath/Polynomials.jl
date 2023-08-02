# Try to keep length based on N,M so no removal of trailing zeros by default
# order is ignored, firstindex is always 0

"""
    ImmutableDensePolynomial{B,T,X,N}

This polynomial type uses an `NTuple{N,T}` to store the coefficients of a polynomial relative to the basis `B` with indeterminate `X`.
For type stability, these polynomials may have trailing zeros. For example, the polynomial `p-p` will have the same size
coefficient tuple as `p`. The `chop` function will trim off trailing zeros, when desired.

Immutable is a bit of a misnomer, as using the `@set!` macro from `Setfield.jl` one can modify elements, as in `@set! p[i] = value`.
"""
struct ImmutableDensePolynomial{B,T,X,N} <: AbstractDenseUnivariatePolynomial{B,T,X}
    coeffs::NTuple{N,T}
    function ImmutableDensePolynomial{B,T,X,N}(cs::NTuple{N,S}) where {B,N,T,X,S}
        new{B,T,Symbol(X),N}(cs)
    end
end

ImmutableDensePolynomial{B,T,X,N}(check::Type{Val{false}}, cs::NTuple{N,T}) where {B,N,T,X} =
    ImmutableDensePolynomial{B,T,X}(cs)

ImmutableDensePolynomial{B,T,X,N}(check::Type{Val{true}}, cs::NTuple{N,T}) where {B,N, T,X} =
    ImmutableDensePolynomial{B,T,X,N}(cs)

# tuple with mis-matched size
function ImmutableDensePolynomial{B,T,X,N}(xs::NTuple{M,S}) where {B,T,S,X,N,M}
    p = ImmutableDensePolynomial{B,S,X,M}(xs)
    convert(ImmutableDensePolynomial{B,T,X,N}, ImmutableDensePolynomial{B,T,X,M}(xs))
end

# constant
function ImmutableDensePolynomial{B,T,X,N}(c::S) where {B,T,X,N,S<:Scalar}
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
ImmutableDensePolynomial{B,T,X}(::Val{false}, xs::NTuple{N,S}) where {B,T,S,X,N} = ImmutableDensePolynomial{B,T,X,N}(convert(NTuple{N,T}, xs))
ImmutableDensePolynomial{B,T,X}(xs::NTuple{N,S}) where {B,T,S,X,N} = ImmutableDensePolynomial{B,T,X,N}(convert(NTuple{N,T}, xs))
ImmutableDensePolynomial{B,T}(xs::NTuple{N,S}, var::SymbolLike=Var(:x)) where {B,T,S,N} = ImmutableDensePolynomial{B,T,Symbol(var),N}(xs)
ImmutableDensePolynomial{B}(xs::NTuple{N,T}, var::SymbolLike=Var(:x)) where {B,T,N} = ImmutableDensePolynomial{B,T,Symbol(var),N}(xs)

# abstract vector. Must eat order
ImmutableDensePolynomial{B,T,X}(::Val{false}, xs::AbstractVector{S}, order::Int=0) where {B,T,X,S} =
    ImmutableDensePolynomial{B,T,X}(xs, order)

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
    ImmutableDensePolynomial{B,T,X,N}(ntuple(i -> T(p[i-1]), Val(N)))
end

function Base.map(fn, p::P, args...)  where {B,T,X, P<:ImmutableDensePolynomial{B,T,X}}
    xs = map(fn, p.coeffs, args...)
    R = eltype(xs)
    return ImmutableDensePolynomial{B,R,X}(xs)
end


Base.copy(p::ImmutableDensePolynomial) = p
Base.similar(p::ImmutableDensePolynomial, args...) = p.coeffs


# not type stable, as N is value dependent
function trim_trailing_zeros(cs::Tuple)
    isempty(cs) && return cs
    !iszero(last(cs)) && return cs
    i = findlast(!iszero, cs)
    i == nothing && return ()
    xs = ntuple(Base.Fix1(getindex,cs), Val(i))
    xs
end

## chop. Also, not type stable
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

# misnamed, should be chop!!
chop!(p::ImmutableDensePolynomial; kwargs...) = chop(p; kwargs...)

# truncate!!; keeps length replacing values with zeros
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


# isapprox helper
function normΔ(q1::ImmutableDensePolynomial{B}, q2::ImmutableDensePolynomial{B}) where {B}
    iszero(q1) && return norm(q2, p)
    iszero(q2) && return norm(q1, p)
    r = abs(zero(q1[end] + q2[end]))
    tot = zero(r)
    for i ∈ 1:maximum(lastindex, (q1,q2))
       @inbounds tot += abs2(q1[i] - q2[i])
    end
    return sqrt(tot)
end


## ---

_zeros(::Type{<:ImmutableDensePolynomial}, z::S, N) where {S} =
    ntuple(_ -> zero(S), Val(N))

minimumexponent(::Type{<:ImmutableDensePolynomial}) =  0

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

# need to call with Setfield as in
# @set! p[i] = value
function Base.setindex(p::ImmutableDensePolynomial{B,T,X,N}, value, i::Int) where {B,T,X,N}
    ps = p.coeffs
    @set! ps[i] = value
    ImmutableDensePolynomial{B,T,X,N}(ps)
end

Base.setindex!(p::ImmutableDensePolynomial, value, i::Int) =
    throw(ArgumentError("Use the `@set!` macro from `Setfield` to mutate coefficients."))


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



## ---
degree(p::ImmutableDensePolynomial{B,T,X,0}) where {B,T,X} = -1
function degree(p::ImmutableDensePolynomial{B,T,X,N}) where {B,T,X,N}
    i = findlast(!iszero, p.coeffs)
    isnothing(i) && return -1
    return i - 1
end

function coeffs(p::P) where {P <: ImmutableDensePolynomial}
    trim_trailing_zeros(p.coeffs)
end

# zero, one
Base.zero(::Type{<:ImmutableDensePolynomial{B,T,X}}) where {B,T,X} =
    ImmutableDensePolynomial{B,T,X,0}(())

function Base.zero(P::Type{ImmutableDensePolynomial{B,T,X,N}}) where {B,T,X,N}
    xs = _zeros(P, zero(T), N)
    ImmutableDensePolynomial{B,T,X,N}(xs)
end

function basis(P::Type{<:ImmutableDensePolynomial{B}}, i::Int) where {B}
    xs = _zeros(P, zero(eltype(P)), i + 1)
    @set! xs[end] = 1
    ImmutableDensePolynomial{B,eltype(P),indeterminate(P)}(xs)
end



## Vector space operations
## vector ops +, -, c*x

## binary
function Base.:+(p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where{B,X,T,S,N,M}
    N < M && return q + p
    _tuple_combine(+, p, q)
end
function Base.:-(p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where{B,X,T,S,N,M}
    N < M && return (-q) + p
    _tuple_combine(-, p, q)
end

# handle +, -; Assume N >= M
_tuple_combine(op, p::ImmutableDensePolynomial{B,T,X,0}, q::ImmutableDensePolynomial{B,S,X,M}) where{B,X,T,S,M} =
    zero(ImmutableDensePolynomial{B,T,X,0})
function _tuple_combine(op, p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where{B,X,T,S,N,M}
    @assert N >= M
    xs = _tuple_combine(op, p.coeffs, q.coeffs)
    R = eltype(xs)
    ImmutableDensePolynomial{B,R,X,N}(xs)
end


# scalar mult faster with 0 default
scalar_mult(p::ImmutableDensePolynomial{B,T,X,0}, c::S) where {B,T,X,S} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})
scalar_mult(c::S, p::ImmutableDensePolynomial{B,T,X,0}) where {B,T,X,S} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})


## ---
# Padded vector combination of two homogeneous tuples assuming N ≥ M
@generated function _tuple_combine(op, p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}

    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(op(p1[$i],p2[$i]))
    end
    for i in (M+1):N
        exprs[i] =:(p1[$i])
    end

    return quote
        Base.@_inline_meta
        #Base.@inline
        tuple($(exprs...))
    end

end

## Static size of product makes generated functions  a good choice
## from https://github.com/tkoolen/StaticUnivariatePolynomials.jl/blob/master/src/monomial_basis.jl
## convolution of two tuples
@generated function fastconv(p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}
    P = M + N - 1
    exprs = Any[nothing for i = 1 : P]
    for i in 1 : N
        for j in 1 : M
            k = i + j - 1
            if isnothing(exprs[k])
                exprs[k] = :(p1[$i] * p2[$j])
            else
                exprs[k] = :(muladd(p1[$i], p2[$j], $(exprs[k])))
            end
        end
    end

    return quote
        Base.@_inline_meta # 1.8 deprecation
        tuple($(exprs...))
    end

end
