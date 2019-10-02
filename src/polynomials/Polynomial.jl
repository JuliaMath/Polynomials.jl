using LinearAlgebra

export Polynomial

const SymbolLike = Union{AbstractString,Char,Symbol}

struct Polynomial{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Polynomial(a::AbstractVector{T}, var::SymbolLike = :x) where {T <: Number}
        # if a == [] we replace it with a = [0]
        if length(a) == 0
            return new{T}(zeros(T, 1), Symbol(var))
        else
          # determine the last nonzero element and truncate a accordingly
            last_nz = findlast(!iszero, a)
            a_last = max(1, last_nz === nothing ? 0 : last_nz)
            new{T}(a[1:a_last], Symbol(var))
        end
    end
end

Polynomial(n::Number, var::SymbolLike = :x) = Polynomial([n], var)
function Polynomial{T}(n::S, var::SymbolLike = :x) where {T,S<:Number}
    U = promote_type(T, S)
    Polynomial{U}([n], var)
end
function Polynomial{T}(x::AbstractVector{S}, var::SymbolLike = :x) where {T <: Number,S <: Number}
    y = convert(Vector{T}, x)
    Polynomial(y, var)
end
Polynomial(x, var::SymbolLike = :x) = Polynomial(collect(x), var)

_domain(p::Polynomial) = (-∞, ∞)
scale_to_domain(P::Type{Polynomial}, x) = x

function (p::Polynomial{T})(x::S) where {T,S}
    R = promote_type(T, S)
    length(p) == 0 && return zero(R)
    y = convert(R, p[end])
    @inbounds for i in (lastindex(p) - 1):-1:0
        y = p[i] + x * y
    end
    return y    
end

function _fromroots(P::Type{Polynomial}, r::AbstractVector{T}, var::SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j = 1:n, i = j:-1:1
        c[i + 1] = c[i + 1] - r[j] * c[i]
    end
    return Polynomial(reverse(c), var)
end


function _vander(P::Type{Polynomial}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    @inbounds for i in 1:n
        A[:, i + 1] = A[:, i] .* x
    end
    return A
end

function _integral(p::Polynomial{T}, k::S) where {T,S}
    n = length(p)
    R = promote_type(typeof(one(T) / 1), S)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return Polynomial(a2, p.var)
end

function _derivative(p::Polynomial{T}, k::Integer) where {T}
    n = length(p)
    a2 = Vector{T}(undef, n - k)
    @inbounds for i in k:n - 1
        a2[i - k + 1] = reduce(*, (i - k + 1):i, init = p[i])
    end
    return Polynomial(a2, p.var)
end

function _companion(p::Polynomial{T}) where T
    d = degree(p)
    comp = diagm(-1 => ones(T, d - 1))
    monics = p.coeffs ./ p.coeffs[end]
    comp[:, end] .= -monics[1:d]
    return comp
end

_mul(p::Polynomial, c) = Polynomial(p.coeffs .* c, p.var)
_div(p::Polynomial, c) = Polynomial(p.coeffs ./ c, p.var)