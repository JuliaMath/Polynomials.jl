## Code from aamini/FastConv.jl for convolutions of small size

using Base.Cartesian


# direct version (do not check if threshold is satisfied)
@generated function fastconv(E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote
        retsize = [size(E)...] + [size(k)...] .- 1
        retsize = tuple(retsize...)
        ret = zeros(T, retsize)
        convn!(ret, E, k)
        return ret
    end
end

# in place helper operation to speedup memory allocations
@generated function convn!(out::Array{T}, E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote
        @inbounds begin
            @nloops $N x E begin
                @nloops $N i k begin
                    (@nref $N out d -> (x_d + i_d - 1)) += (@nref $N E x) * (@nref $N k i)
                end
            end
        end
        return out
    end
end


## Code from Julia 1.4   (https://github.com/JuliaLang/julia/blob/master/base/math.jl#L101 on 4/8/20)
## cf. https://github.com/JuliaLang/julia/pull/32753
## Slight modification when `x` is a matrix
## Remove once dependencies for Julia 1.0.0 are dropped
module EvalPoly
using LinearAlgebra
function evalpoly(x::S, p::Tuple) where {S}
    p == () && return zero(S)
    if @generated
        N = length(p.parameters)
        ex = :(p[end]*_one(S))
        for i in N-1:-1:1
            ex = :(_muladd($ex, x, p[$i]))
        end
        ex
    else
        _evalpoly(x, p)
    end
end

evalpoly(x, p::AbstractVector) = _evalpoly(x, p)

function _evalpoly(x::S, p) where {S}
    N = length(p)
    ex = p[end]*_one(x)
    for i in N-1:-1:1
        ex = _muladd(ex, x, p[i])
    end
    ex
end

function evalpoly(z::Complex, p::Tuple)
    if @generated
        N = length(p.parameters)
        a = :(p[end]*_one(z))
        b = :(p[end-1]*_one(z))
        as = []
        for i in N-2:-1:1
            ai = Symbol("a", i)
            push!(as, :($ai = $a))
            a = :(muladd(r, $ai, $b))
            b = :(p[$i] - s * $ai)
        end
        ai = :a0
        push!(as, :($ai = $a))
        C = Expr(:block,
                 :(x = real(z)),
                 :(y = imag(z)),
                 :(r = x + x),
                 :(s = muladd(x, x, y*y)),
                 as...,
                 :(muladd($ai, z, $b)))
    else
        _evalpoly(z, p)
    end
end
evalpoly(z::Complex, p::Tuple{<:Any}) = p[1]*_one(z)


evalpoly(z::Complex, p::AbstractVector) = _evalpoly(z, p)

function _evalpoly(z::Complex, p)
    length(p) == 1 && return p[1]*_one(z)
    N = length(p)
    a = p[end]*_one(z)
    b = p[end-1]*_one(z)

    x = real(z)
    y = imag(z)
    r = 2x
    s = muladd(x, x, y*y)
    for i in N-2:-1:1
        ai = a
        a = muladd(r, ai, b)
        b = p[i] - s * ai
    end
    ai = a
    muladd(ai, z, b)
end

## modify muladd, as needed
_muladd(a,b,c) = muladd(a,b,c)
_muladd(a, b::Vector, c) = a.*b .+ c
_muladd(a, b::Matrix, c) = (a*I)*b + c*I

# try to get y = P(c::T)(x::S) = P{T}(c)(x::S) to
# have y == one(T)*one(S)*x
_one(P::Type{<:Matrix}) = one(eltype(P))*I
_one(x::Matrix) = one(eltype(x))*I
_one(x) = one(x)
end

## get type of parametric composite type without type parameters
## this is needed when the underlying type changes, e.g. with integration
## where T=Int might become T=Float64
##
## trick from [ConstructionBase.jl](https://github.com/JuliaObjects/ConstructionBase.jl/blob/b5686b755bd3bee29b181b3cb18fe2effa0f10a2/src/ConstructionBase.jl#L25)
## as noted in https://discourse.julialang.org/t/get-new-type-with-different-parameter/37253/4
##
#@generated function constructorof(::Type{T}) where T
#    getfield(parentmodule(T), nameof(T))
#end

# https://discourse.julialang.org/t/how-do-a-i-get-a-type-stripped-of-parameters/73465/11
constructorof(::Type{T}) where T = Base.typename(T).wrapper


# Define our own minimal Interval type, inspired by Intervals.jl.
# We vendor it in to avoid adding the heavy Intervals.jl dependency.
abstract type Bound end
abstract type Bounded <: Bound end
struct Closed <: Bounded end
struct Open <: Bounded end
struct Unbounded <: Bound end

"""
    Interval{T, L <: Bound, R <: Bound}

Very bare bones Interval following `Intervals.jl` assuming `T<:Real`.
"""

struct Interval{T, L <: Bound, R <: Bound}
    first::T
    last::T
    function Interval{T,L,R}(f::T, l::T) where {T, L <: Bound, R <: Bound}
        f < l && return new{T,L,R}(f, l)
        throw(ArgumentError("first not less than last"))
    end
end
function Base.show(io::IO, I::Interval{T,L,R}) where {T,L,R}
    l,r = extrema(I)
    print(io, L == Closed ? "[" : "(")
    print(io, l, ", ", r)
    print(io, R == Closed ? "]" : ")")
end

# type is "[], (), (] or [)"
const _interval_types =
    Dict("[]" => (Closed, Closed), "()" => (Open,Open),
         "[)" => (Closed, Open),   "(]" => (Open, Closed))
"""
    Polynomials.Interval(f, l, typ="[]")

Constructor for return type of `Polynomials.domain`. Default is a closed interval, unless values are infinite.
"""
function Interval(f,l, typ="[]")
    𝐟,𝐥 = promote(f,l)
    L,R = _interval_types[typ]
    𝐋 = isinf(𝐟) ? Unbounded : L
    𝐑 = isinf(𝐥) ? Unbounded : R
    Interval{typeof(𝐟),𝐋,𝐑}(f,l)
end

intervaltype(I::Interval{T,L,R}) where {T,L,R} = (L,R)

Base.first(I::Interval) = I.first
Base.last(I::Interval) = I.last
Base.extrema(I::Interval) = (first(I), last(I))

function Base.in(x, I::Interval{T,L,R}) where {T, L, R}
    a, b = extrema(I)
    (L == Open ? a < x : a <= x) && (R == Open ? x < b : x <= b)
end

Base.isopen(I::Interval{T,L,R}) where {T,L,R} = (L != Closed && R != Closed)
isclosed(I::Interval{T,L,R}) where {T,L,R} = (L == Closed && R == Closed)
isbounded(I::Interval{T,L,R}) where {T,L,R} = (L != Unbounded && R != Unbounded)
