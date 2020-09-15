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


## Code from https://github.com/JuliaMatrices/SpecialMatrices.jl/blob/master/src/vandermonde.jl
## Avoids dependency for fitting Vandermonde(xs) \ ys
## (also Algorithm 2 of https://www.maths.manchester.ac.uk/~higham/narep/narep108.pdf)
"""
    dvand!(a, b) -> b
Solves system ``A*x = b`` in-place.
``A`` is Vandermonde matrix ``A_{ij} = a_i^{j-1}``.
Algorithm by Bjorck & Pereyra,
Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903,
https://doi.org/10.2307/2004623
"""
function dvand!(alpha, B)
    n = length(alpha)
    if n != size(B,1)
        throw(DimensionMismatch("matrix has dimensions ($n,$n) but right hand side has $(size(B,1)) rows"))
    end
    nrhs = size(B,2)
    @inbounds begin
        for j=1:nrhs
            for k=1:n-1
                for i=n:-1:k+1
                    B[i,j] = (B[i,j]-B[i-1,j])/(alpha[i]-alpha[i-k])
                end
            end
            for k=n-1:-1:1
                for i=k:n-1
                    B[i,j] = B[i,j]-alpha[k]*B[i+1,j]
                end
            end
        end
    end
end

"""
    pvand!(a, b) -> b
Solves system ``A^T*x = b`` in-place.
``A^T`` is transpose of Vandermonde matrix ``A_{ij} = a_i^{j-1}``.
Algorithm by Bjorck & Pereyra,
Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903,
https://doi.org/10.2307/2004623
"""
function pvand!(alpha, B)
    n = length(alpha);
    if n != size(B,1)
        throw(DimensionMismatch("matrix has dimensions ($n,$n) but right hand side has $(size(B,1)) rows"))
    end
    nrhs = size(B,2)
    @inbounds begin
        for j=1:nrhs
            for k=1:n-1
                for i=n:-1:k+1
                    B[i,j] = B[i,j]-alpha[k]*B[i-1,j]
                end
            end
            for k=n-1:-1:1
                for i=k+1:n
                    B[i,j] = B[i,j]/(alpha[i]-alpha[i-k])
                end
                for i=k:n-1
                    B[i,j] = B[i,j]-B[i+1,j]
                end
            end
        end
    end
end


## Code from Julia 1.4   (https://github.com/JuliaLang/julia/blob/master/base/math.jl#L101 on 4/8/20)
## cf. https://github.com/JuliaLang/julia/pull/32753
## Slight modification when `x` is a matrix
## Remove once dependencies for Julia 1.0.0 are dropped
function evalpoly(x::S, p::Tuple) where {S}
    p == () && return zero(S)
    if @generated
        N = length(p.parameters)
        ex = :(p[end]*_one(S))
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
    ex = p[end]*_one(x)
    for i in N-1:-1:1
        ex = _muladd(x, ex, p[i])
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
_muladd(a::Vector, b, c) = a.*b .+ c
_muladd(a::Matrix, b, c) = a*(b*I) + c*I

# try to get y = P(c::T)(x::S) = P{T}(c)(x::S) to
# have y == one(T)*one(S)*x
_one(P::Type{<:Matrix}) = one(eltype(P))*I
_one(x::Matrix) = one(eltype(x))*I
_one(x) = one(x)

## get type of parametric composite type without type parameters
## this is needed when the underlying type changes, e.g. with integration
## where T=Int might become T=Float64
##
## trick from [ConstructionBase.jl](https://github.com/JuliaObjects/ConstructionBase.jl/blob/b5686b755bd3bee29b181b3cb18fe2effa0f10a2/src/ConstructionBase.jl#L25)
## as noted in https://discourse.julialang.org/t/get-new-type-with-different-parameter/37253/4
##
@generated function constructorof(::Type{T}) where T
    getfield(parentmodule(T), nameof(T))
end

