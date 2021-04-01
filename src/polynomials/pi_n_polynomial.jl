"""
    ΠₙPolynomial{T,X,N}(coeffs::Vector{T})

Construct a polynomial in `Πₙ`, the collection of polynomials of degree `n` or less using a vector of length `N+1`.

* Unlike other polynomial types, this does not copy the coefficients on construction
* Unlike other polynomial types, this type broadcasts like a vector for in-place vector operations (scalare multiplication, polynomial addition/subtraction of the same size)

"""
struct ΠₙPolynomial{T,X,N} <: Polynomials.StandardBasisPolynomial{T, X}
    coeffs::Vector{T}
    function ΠₙPolynomial{T, X, N}(coeffs::AbstractVector{T}) where {T, X, N}
        n = length(coeffs)
        n > N+1 && throw(ArgumentError("Too many coefficient"))
        n < N+1 && append!(coeffs, zeros(T, N+1-n))
        new{T,X, N}(coeffs) 
    end
    function ΠₙPolynomial{T, X}(coeffs::AbstractVector{T}) where {T, X}
        N = length(coeffs) - 1
        new{T,X,N}(coeffs) # NO CHECK on trailing zeros
    end
    function ΠₙPolynomial(coeffs::AbstractVector{T}, var=:x) where {T}
        N = length(coeffs) - 1
        X = Symbol(var)
        new{T,X,N}(coeffs)
    end
end



Polynomials.@register ΠₙPolynomial

Base.broadcastable(p::ΠₙPolynomial) = p.coeffs;
Base.ndims(::Type{<:ΠₙPolynomial}) = 1
Base.copyto!(p::ΠₙPolynomial, x) = copyto!(p.coeffs, x);

Πₙ(::Type{P}) where {T,X,N, P<:ΠₙPolynomial{T,X,N}} = N # get N
Πₙ(p::ΠₙPolynomial{T,X,N}) where {T,X,N} = N
function Polynomials.degree(p::ΠₙPolynomial)
    i = findlast(!iszero, p.coeffs)
    i == nothing && return -1
    i - 1
end

Polynomials.zero(::Type{P}) where {T,X,N,P <: ΠₙPolynomial{T,X,N}} = P(zeros(T,N+1))
Polynomials.zero(p::P) where {T,X,N,P <: ΠₙPolynomial{T,X,N}} = zero(P)
function Polynomials.one(::Type{P}) where {T,X,N,P <: ΠₙPolynomial{T,X,N}}
    cs = zeros(T, N+1)
    cs[1] = one(T)
    P(cs)
end
Polynomials.one(p::P) where {T,X,N,P} = one(P)

function Polynomials.variable(::Type{P}) where {T,X,N,P <: ΠₙPolynomial{T,X,N}}
    cs = zeros(T, N+1)
    cs[2] = one(T)
    P(cs)
end
Polynomials.variable(p::P) where {T,X,N,P <: ΠₙPolynomial{T,X,N}} = variable(P)


# pre-allocated multiplication
function LinearAlgebra.mul!(pq, p::ΠₙPolynomial{T,X}, q) where {T,X}
    m,n = degree(p), degree(q)
    pq.coeffs .= zero(T)
    for i ∈ 0:m
        for j ∈ 0:n
            k = i + j
            pq.coeffs[1+k] += p.coeffs[1+i] * q.coeffs[1+j]
        end
    end
    nothing
end

# also just p .*= -1
function Base.:-(p::ΠₙPolynomial{T,X}) where {T,X}
    for i ∈ eachindex(p.coeffs)
        p.coeffs[i] *= -1
    end
    p
end



