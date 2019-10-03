using LinearAlgebra

export Polynomial

struct Polynomial{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Polynomial{T}(coeffs::Vector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

Polynomial(a::AbstractVector{T}, var::SymbolLike = :x) where {T} = Polynomial{T}(a, Symbol(var))

Polynomial(n::Number, var = :x) = Polynomial([n], var)
Polynomial{T}(n::S, var = :x) where {T,S <: Number} = Polynomial(T(n), var)
Polynomial{T}(x::AbstractVector{S}, var = :x) where {T,S <: Number} = Polynomial(T.(x), var)
Polynomial(x, var = :x) = Polynomial(collect(x), var)

@register Polynomial


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

function _integral(p::Polynomial{T}, k::S) where {T,S <: Number}
    n = length(p)
    R = promote_type(eltype(one(T) / 1), S)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return Polynomial(a2, p.var)
end

function _derivative(p::Polynomial{T}, order::Integer) where {T}
    n = length(p)
    a2 = Vector{T}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return Polynomial(a2, p.var)
end

function _companion(p::Polynomial{T}) where T
    d = degree(p)
    R = eltype(one(T) / p.coeffs[end])
    comp = diagm(-1 => ones(R, d - 1))
    monics = p.coeffs ./ p.coeffs[end]
    comp[:, end] .= -monics[1:d]
    return comp
end
