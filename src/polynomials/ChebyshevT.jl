export ChebyshevT

"""
    ChebyshevT{<:Number}(coeffs::AbstractVector, var=:x)

Chebyshev polynomial of the first kind
"""
struct ChebyshevT{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function ChebyshevT{T}(coeffs::AbstractVector{T}, var::Symbol) where {T<:Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

@register ChebyshevT

domain(::Type{<:ChebyshevT}) = Interval(-1, 1)
scale_to_domain(::Type{<:ChebyshevT}, x) = @. x * 2 / sum(x) - 1

"""
    (::ChebyshevT)(x)

Evaluate the Chebyshev polynomial at `x`. If `x` is outside of the domain of [-1, 1], an error will be thrown. THe evaluation uses Clenshaw recursion, or synthetic division.
"""
function (ch::ChebyshevT{T})(x::S) where {T,S}
    R = promote_type(T, S)
    length(ch) == 0 && return zero(R)
    if length(ch) < 3
        c0 = ch[0]
        c1 = ch[1]
    else
        c0 = ch[0]
        c1 = ch[1]
        for i in 2:length(ch) - 1
            c0, c1 = ch[i] - c1, c0 + c1 * 2x
        end
    end
    return R(c0 + c1 * x)
end
#####
function fromroots(P::Type{<:ChebyshevT}, r::AbstractVector{T}; var::SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j = 1:n, i = j:-1:1
        c[i + 1] = c[i + 1] - r[j] * c[i]
    end
    return ChebyshevT(reverse(c), var)
end

function vander(P::Type{<:ChebyshevT}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    @inbounds for i in 1:n
        A[:, i + 1] = A[:, i] .* x
    end
    return A
end

function integral(p::ChebyshevT{T}, k::S) where {T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(k)
        return ChebyshevT([NaN])
    end
    n = length(p)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return ChebyshevT(a2, p.var)
end

function derivative(p::ChebyshevT{T}, order::Integer) where {T}
    order < 0 && error("Order of derivative must be non-negative")
    order == 0 && return p
    hasnan(p) && return ChebyshevT(T[NaN], p.var)
    order > length(p) && return zero(ChebyshevT{T})

    n = length(p)
    a2 = Vector{T}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return ChebyshevT(a2, p.var)
end

function companion(p::ChebyshevT{T}) where T
    d = length(p) - 1
    d < 1 && error("Series must have degree greater than 1")
    d == 1 && return diagm(0 => [-p[0] / p[1]])

    R = eltype(one(T) / p.coeffs[end])
    comp = diagm(-1 => ones(R, d - 1))
    monics = p.coeffs ./ p.coeffs[end]
    comp[:, end] .= -monics[1:d]
    return comp
end

function Base.:+(p1::ChebyshevT, p2::ChebyshevT)
    p1.var != p2.var && error("Polynomials must have same variable")
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n]
    return ChebyshevT(c, p1.var)
end

function Base.:*(p1::ChebyshevT{T}, p2::ChebyshevT{S}) where {T,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    n = length(p1) - 1
    m = length(p2) - 1
    R = promote_type(T, S)
    c = zeros(R, m + n + 1)
    for i = 0:n, j = 0:m
        c[i + j + 1] += p1[i] * p2[j]
    end
    return ChebyshevT(c, p1.var)
end

function Base.divrem(num::ChebyshevT{T}, den::ChebyshevT{S}) where {T,S}
    num.var != den.var && error("Polynomials must have same variable")
    n = length(num) - 1
    m = length(den) - 1
    if m == 0 && den[0] ≈ 0 throw(DivideError()) end
    R = typeof(one(T) / one(S))
    P = ChebyshevT{R}
    deg = n - m + 1
    if deg ≤ 0 
        return zero(P), convert(P, num)
    end
    q_coeff = zeros(R, deg)
    r_coeff = R.(num[0:n])
    @inbounds for i in n:-1:m
        q = r_coeff[i + 1] / den[m]
        q_coeff[i - m + 1] = q
        @inbounds for j in 0:m
            elem = den[j] * q
            r_coeff[i - m + j + 1] -= elem
        end
    end
    return P(q_coeff, num.var), P(r_coeff, num.var)
end

function Base.gcd(a::ChebyshevT{T}, b::ChebyshevT{S}) where {T,S}
  U       = typeof(one(T) / one(S))
  r₀ = convert(ChebyshevT{U}, a)
  r₁ = truncate!(convert(ChebyshevT{U}, b))
  iter    = 1
  itermax = length(b)

  while r₁ ≉ zero(r₁) && iter ≤ itermax   # just to avoid unnecessary recursion
    _, rtemp  = divrem(r₀, r₁)
    r₀        = r₁
    r₁        = truncate(rtemp)
    iter      += 1
  end
  return r₀
end

function printpoly(io::IO, p::ChebyshevT, mimetype=MIME"text/plain"(); descending_powers=false, offset::Int=0)
    chopped = chop(p)
    print(io, coeffs(chopped))
    return nothing
end