export Polynomial

const SymbolLike = Union{AbstractString,Char,Symbol}

struct Polynomial{T <: Number} <: AbstractPolynomial
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

domain(p::Polynomial) = (-∞, ∞)

function _fromroots(P::Type{Polynomial}, r::AbstractVector{T}, var::SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j = 1:n, i = j:-1:1
        c[i + 1] = c[i + 1] - r[j] * c[i]
    end
    return Polynomial(reverse(c), var)
end
