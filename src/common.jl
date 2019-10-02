using LinearAlgebra
import Base: ==
import InteractiveUtils: subtypes

export AbstractPolynomial,
       fromroots,
       truncate!,
       chop!,
       coeffs,
       degree,
       hasnan,
       roots,
       companion,
       vander,
       fit,
       integrate,
       integral,
       derivative


const SymbolLike = Union{AbstractString,Char,Symbol}
abstract type AbstractPolynomial{T <: Number} end

"""
    fromroots(::AbstractVector{<:Number}, var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractVector{<:Number}, var=:x)

Construct a polynomial of the given type given the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest
julia> r = [3, 2]; # (x - 3)(x - 2)

julia> fromroots(r)
Polynomial(x^2 - 5x + 6)
```
"""
fromroots(P::Type{<:AbstractPolynomial}, r::AbstractVector{T}, var::SymbolLike = :x) where {T <: Number} = _fromroots(P, r, var)
fromroots(r::AbstractVector{<:Number}, var::SymbolLike = :x) = fromroots(Polynomial, r, var)
fromroots(r, var::SymbolLike = :x) = fromroots(collect(r), var)

"""
    fromroots(::AbstractMatrix{<:Number}, var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractMatrix{<:Number}, var=:x)

Construct a polynomial of the given type using the eigenvalues of the given matrix as the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest
julia> A = [1 2; 3 4]; # (x - 5.37228)(x + 0.37228)

julia> fromroots(A)
Polynomial(-1.9999999999999998 - 5.0⋅x + 1.0⋅x^2)
```
"""
fromroots(P::Type{<:AbstractPolynomial}, A::AbstractMatrix{T}, var::SymbolLike = :x) where {T <: Number} = fromroots(P, eigvals(A), var)
fromroots(A::AbstractMatrix{T}, var::SymbolLike = :x) where {T <: Number} = fromroots(Polynomial, eigvals(A), var)

"""
    fit(x, y, deg=length(x) - 1)
    fit(::Type{<:AbstractPolynomial}, x, y, deg=length(x) - 1)

Fit the given data as a polynomial type with the given degree. The default polynomial type is [`Polynomial`](@ref).
"""
fit(x, y, deg = length(x) - 1) = fit(Polynomial, collect(x), collect(y), deg)
fit(P::Type{<:AbstractPolynomial}, x, y, deg = length(x) - 1) = fit(P, collect(x), collect(y), deg)
function fit(P::Type{<:AbstractPolynomial}, x::AbstractVector, y::AbstractVector, deg::Integer = length(x) - 1)
    x = scale_to_domain(P, x)
    vand = _vander(P, x, deg)
    coeffs = pinv(vand) * y
    return P(coeffs)
end


"""
    roots(::AbstractPolynomial)

Returns the roots of the given polynomial. This is calculated via the eigenvalues of the companion matrix.
"""
function roots(p::AbstractPolynomial{T}) where {T <: Number}
    d = degree(p)
    d < 1 && error("Series must have degree greater than 1")
    d == 1 && return diagm([-p[0] / p[1]])

    chopped_trimmed = chop(truncate(p))
    n_trail = length(p) - length(chopped_trimmed)
    comp = _companion(chopped_trimmed)
    L = eigvals(rot180(comp))
    append!(L, zeros(n_trail))
    return sort!(L, rev = true)
end

"""
    companion(::AbstractPolynomial)

Return the companion matrix for the given polynomial.

# References
[Companion Matrix](https://en.wikipedia.org/wiki/Companion_matrix)
"""
function companion(p::AbstractPolynomial)
    d = degree(p)
    d < 1 && error("Series must have degree greater than 1")
    d == 1 && return diagm([-p[0] / p[1]])
    _companion(p)
end

"""
    vander(::AbstractPolynomial, x::AbstractVector, deg::Integer)

Calculate the psuedo-Vandermonde matrix of the given polynomial type with the given degree.

# References
[Vandermonde Matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix)
"""
vander(P::Type{<:AbstractPolynomial}, x::AbstractVector, deg::Integer) = _vander(P, x, deg)
vander(P::Type{<:AbstractPolynomial}, x, deg::Integer) = vander(P, collect(x), deg)

"""
    integral(::AbstractPolynomial, k=0)

Returns a polynomial that is the integral of the given polynomial with constant term `k` added. 
"""
integral(p::AbstractPolynomial, k::Number = 0) = _integral(p, k)

"""
    integrate(::AbstractPolynomial, a, b)

Compute the definite integral of the given polynomial from `a` to `b`. Will throw an error if either `a` or `b` are out of the polynomial's domain.
"""
function integrate(p::AbstractPolynomial, a::Number, b::Number)
    P = integral(p)
    return P(b) - P(a)
end

"""
    derivative(::AbstractPolynomial, k = 1)

Returns a polynomail that is the `k`th derivative of the given polynomial. `k` must be non-negative.
"""
function derivative(p::P, k = 1) where {P <: AbstractPolynomial}
    k < 0 && error("Order of derivative must be non-negative")
    k == 0 && return p
    k > length(p) && return zero(P)
    _derivative(p, k)
end

"""
    truncate!(::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

In-place version of [`truncate`](@ref)
"""
function truncate!(p::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0) where {T}
    max_coeff = maximum(abs, p.coeffs)
    thresh = max_coeff * rtol + atol
    map!(c->abs(c) <= thresh ? zero(T) : c, p.coeffs, p.coeffs)
    return chop!(p, rtol = rtol, atol = atol)
end

"""
    truncate(::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

Rounds off coefficients close to zero, as determined by `rtol` and `atol`, and then chops any leading zeros. Returns a new polynomial.
"""
function Base.truncate(p::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0) where {T}
    truncate!(deepcopy(p), rtol = rtol, atol = atol)
end

"""
    chop!(::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0))

In-place version of [`chop`](@ref)
"""
function chop!(p::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0) where {T}
    for i in lastindex(p):-1:0
        val = p[i]
        if !isapprox(val, zero(T); rtol = rtol, atol = atol)
            resize!(p.coeffs, i + 1)
            return p
        end
    end 
    resize!(p.coeffs, 1)
    return p
end

"""
    chop(::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0))

Removes any leading coefficients that are approximately 0 (using `rtol` and `atol`). Returns a polynomial whose degree will guaranteed to be equal to or less than the given polynomial's. 
"""
function Base.chop(p::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0) where {T}
    chop!(deepcopy(p), rtol = rtol, atol = atol)
end

#=
Linear Algebra
=#
LinearAlgebra.norm(q::AbstractPolynomial, p::Real = 2) = norm(coeffs(q), p)
LinearAlgebra.conj(p::P) where {P <: AbstractPolynomial} = P(conj(coeffs(p)))

#=
Conversions
=#
Base.convert(::Type{P}, p::P) where {T,P <: AbstractPolynomial{T}} = p
for i in subtypes(AbstractPolynomial)
    @eval begin
        Base.convert(::Type{$i{T}}, p::$i) where {T} = $i(T.(p.coeffs), p.var)
        Base.convert(::Type{$i{T}}, x, var::SymbolLike = :x) where {T} = $i(T.(x), var)
        Base.promote_rule(::Type{$i{T}}, ::Type{$i{S}}) where {T,S} = $i{promote_type(T, S)}
end
end
Base.promote_rule(::Type{<:AbstractPolynomial{T}}, ::Type{<:AbstractPolynomial{S}}) where {T,S} = Polynomial{promote_type(T, S)}

#=
Inspection
=#
Base.length(p::AbstractPolynomial) = length(p.coeffs)
Base.size(p::AbstractPolynomial) = size(p.coeffs)
Base.size(p::AbstractPolynomial, i::Integer) = size(p.coeffs, i)
Base.eltype(p::AbstractPolynomial{T}) where {T} = T
Base.eltype(::Type{<:AbstractPolynomial{T}}) where {T} = T
Base.iszero(p::AbstractPolynomial) = length(p) == 0

"""
    coeffs(::AbstractPolynomial)

Return the coefficient vector `[a_0, a_1, ..., a_n]` of a polynomial.
"""
coeffs(p::AbstractPolynomial) = p.coeffs

"""
    degree(::AbstractPolynomial)

Return the degree of the polynomial, i.e. the highest exponent in the polynomial that
has a nonzero coefficient. The degree of the zero polynomial is defined to be -1.
"""
degree(p::AbstractPolynomial) = iszero(p) ? -1 : length(p) - 1
hasnan(p::AbstractPolynomial) = any(isnan.(p.coeffs))

#=
indexing
=#
Base.firstindex(p::AbstractPolynomial) = 0
Base.lastindex(p::AbstractPolynomial) = degree(p)
Base.getindex(p::AbstractPolynomial{T}, idx::Int) where {T} = (idx ≥ length(p) ? zero(T) : p.coeffs[idx + 1])
Base.getindex(p::AbstractPolynomial, indices::AbstractVector{Int}) = map(idx->p[idx], indices)
Base.getindex(p::AbstractPolynomial, ::Colon) = p[0:length(p) - 1]
Base.broadcastable(p::AbstractPolynomial) = Ref(p)
# Base.collect(p::P) where {P <: AbstractPolynomial} = collect(P, p)
# Base.iterate(p::AbstractPolynomial) = (p[0] * one(typeof(p)), 1)
# Base.iterate(p::AbstractPolynomial, state) = state <= degree(p) ? (p[state] * variable(p)^(state), state + 1) : nothing

function Base.setindex!(p::AbstractPolynomial, value, idx::Int)
    n = length(p)
    if n ≤ idx
        resize!(p.coeffs, idx + 1)
        p.coeffs[n + 1:idx] .= 0
    end
    p.coeffs[idx + 1] = value
    return p
end
function Base.setindex!(p::AbstractPolynomial, values, indices::AbstractVector{Int})
    setindex!.(p, values, indices)
    return p
end
Base.setindex!(p::AbstractPolynomial, values, ::Colon) = setindex!(p, values, 0:degree(p))
Base.eachindex(p::AbstractPolynomial) = 0:degree(p)

#=
identity
=#
Base.copy(p::P) where {P <: AbstractPolynomial} = P(copy(p.coeffs), p.var)
Base.hash(p::AbstractPolynomial, h::UInt) = hash(p.var, hash(p.coeffs, h))
Base.zero(::Type{P}) where {P <: AbstractPolynomial} = P(zeros(1))
Base.one(::Type{P}) where {P <: AbstractPolynomial} = P(ones(1))

#=
arithmetic
=#
Base.:-(p::P) where {P <: AbstractPolynomial} = P(-p.coeffs, p.var) 
Base.:+(c::Number, p::AbstractPolynomial) = +(p, c)
Base.:-(p::AbstractPolynomial, c::Number) = +(p, -c)
Base.:*(c::Number, p::AbstractPolynomial) = *(p, c)


Base.:+(p::AbstractPolynomial, c) = _add(p, c)
Base.:*(p::AbstractPolynomial, c) = _mul(p, c)
Base.:/(p::AbstractPolynomial, c) = _div(p, c)
Base.:-(p1::AbstractPolynomial, p2::AbstractPolynomial) = +(p1, -p2)

function Base.:+(p1::P, p2::O) where {P <: AbstractPolynomial,O <: AbstractPolynomial}
    p1.var != p2.var && error("Polynomials must have same variable")
    T = promote_type(P, O)
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n]
    return T(c, p1.var)
end

function Base.:*(p1::P, p2::O) where {P <: AbstractPolynomial,O <: AbstractPolynomial}
    p1.var != p2.var && error("Polynomials must have same variable")
    n = degree(p1)
    m = degree(p2)
    T = promote_type(P, O)
    c = zeros(m + n + 1)
    for i = 0:n, j = 0:m
        c[i + j + 1] += p1[i] * p2[j]
    end
    return T(c, p1.var)
end

Base.:^(p::AbstractPolynomial, n::Integer) = Base.power_by_squaring(p, n)

#=
Comparisons
=#
Base.isequal(p1::P, p2::P) where {P <: AbstractPolynomial} = hash(p1) == hash(p2)
==(p1::P, p2::P) where {P <: AbstractPolynomial} = p1.coeffs == p2.coeffs && p1.var == p2.var
==(p1::AbstractPolynomial, n::Number)  = p1.coeffs == [n]
==(n::Number, p1::AbstractPolynomial)  = p1 == n

function Base.isapprox(p1::AbstractPolynomial{T}, p2::AbstractPolynomial{S}; rtol::Real = (Base.rtoldefault(T, S, 0)), atol::Real = 0) where {T,S} 
    isapprox(p1.coeffs, p2.coeffs, rtol = rtol, atol = atol)
end

function Base.isapprox(p1::AbstractPolynomial{T}, n::S; rtol::Real = (Base.rtoldefault(T, S, 0)), atol::Real = 0) where {T,S}
    isapprox(p1.coeffs, [n], rtol = rtol, atol = atol)
end

Base.isapprox(n::S, p1::AbstractPolynomial{T}; rtol::Real = (Base.rtoldefault(T, S, 0)), atol::Real = 0) where {T,S} = isapprox(p1, n, rtol = rtol, atol = atol)
