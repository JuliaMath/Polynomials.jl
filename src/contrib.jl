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
## Remove once dependencies for Julia 1.0.0
function evalpoly(x::S, p::Tuple) where {S}
    if @generated
        N = length(p.parameters)
        ex = :(p[end])
        for i in N-1:-1:1
            ex = :(_muladd(x, $ex, p[$i]))
        end
        ex
    else
        _evalpoly(x, p)
    end
end

evalpoly(x, p::AbstractVector) = _evalpoly(x, p)

function _evalpoly(x::S, p) where {S}
    N = length(p)
    ex = p[end]
    for i in N-1:-1:1
        ex = _muladd(x, ex, p[i])
    end
    ex
end

function evalpoly(z::Complex, p::Tuple)
    if @generated
        N = length(p.parameters)
        a = :(p[end])
        b = :(p[end-1])
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
evalpoly(z::Complex, p::Tuple{<:Any}) = p[1]


evalpoly(z::Complex, p::AbstractVector) = _evalpoly(z, p)

function _evalpoly(z::Complex, p)
    length(p) == 1 && return p[1]
    N = length(p)
    a = p[end]
    b = p[end-1]

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

## modify muladd for matrices
_muladd(a,b,c) = muladd(a,b,c)
_muladd(a::Matrix, b, c) = a*(b*I) + c*I

