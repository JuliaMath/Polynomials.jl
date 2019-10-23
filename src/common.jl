using LinearAlgebra

export fromroots,
       truncate!,
       chop!,
       coeffs,
       degree,
       domain,
       scale_to_domain,
       order,
       hasnan,
       roots,
       companion,
       vander,
       fit,
       integrate,
       integral,
       derivative,
       variable

"""
    fromroots(::AbstractVector{<:Number}; var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractVector{<:Number}; var=:x)

Construct a polynomial of the given type given the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest
julia> r = [3, 2]; # (x - 3)(x - 2)

julia> fromroots(r)
Polynomial(6 - 5*x + x^2)
```
"""
fromroots(P::Type{<:AbstractPolynomial}, r::AbstractVector; var::SymbolLike = :x)
fromroots(r::AbstractVector{<:Number}; var::SymbolLike = :x) = fromroots(Polynomial, r, var = var)
fromroots(r; var::SymbolLike = :x) = fromroots(collect(r), var = var)

"""
    fromroots(::AbstractMatrix{<:Number}; var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractMatrix{<:Number}; var=:x)

Construct a polynomial of the given type using the eigenvalues of the given matrix as the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest
julia> A = [1 2; 3 4]; # (x - 5.37228)(x + 0.37228)

julia> fromroots(A)
Polynomial(-1.9999999999999998 - 5.0*x + 1.0*x^2)
```
"""
fromroots(P::Type{<:AbstractPolynomial}, A::AbstractMatrix{T}; var::SymbolLike = :x) where {T <: Number} = fromroots(P, eigvals(A), var = var)
fromroots(A::AbstractMatrix{T}; var::SymbolLike = :x) where {T <: Number} = fromroots(Polynomial, eigvals(A), var = var)

"""
    fit(x, y; [weights], deg=length(x) - 1, var=:x)
    fit(::Type{<:AbstractPolynomial}, x, y; [weights], deg=length(x)-1, var=:x)

Fit the given data as a polynomial type with the given degree. Uses linear least squares. When weights are given, as either a `Number`, `Vector` or `Matrix`, will use weighted linear least squares. The default polynomial type is [`Polynomial`](@ref).
"""
function fit(P::Type{<:AbstractPolynomial}, 
            x::AbstractVector{T}, 
            y::AbstractVector{T}; 
            weights = nothing, deg::Integer = length(x) - 1, var = :x) where {T <: Number}
    vand = vander(P, x, deg)
    if weights !== nothing
        coeffs = _wlstsq(vand, y, weights)
    else
        coeffs = pinv(vand) * y
    end
    return P(T.(coeffs), var)
end

fit(x, y; weights = nothing, deg = length(x) - 1, var = :x) = fit(Polynomial{eltype(y)}, collect(x), collect(y); weights = weights, deg = deg, var = var)

fit(P::Type{<:AbstractPolynomial}, x, y; weights = nothing, deg = length(x) - 1, var = :x) = fit(P, collect(x), collect(y); weights = weights, deg = deg, var = var)

function fit(P::Type{<:AbstractPolynomial}, 
    x::AbstractVector{T}, 
    y::AbstractVector{S}; 
    weights = nothing, deg = length(x) - 1, var = :x) where {T,S}
    x, y = promote(x, y)
    fit(P, x, y; weights = weights, deg = deg, var = var)
end

# Weighted linear least squares
_wlstsq(vand, y, W::Number) = _wlstsq(vand, y, fill!(similar(y), W))
_wlstsq(vand, y, W::AbstractVector) = _wlstsq(vand, y, diagm(0 => W))
_wlstsq(vand, y, W::AbstractMatrix) = (vand' * W * vand) \ (vand' * W * y)

"""
    roots(::AbstractPolynomial)

Returns the roots of the given polynomial. This is calculated via the eigenvalues of the companion matrix.
"""
function roots(p::AbstractPolynomial{T}) where {T <: Number}
    d = length(p) - 1
    if d < 1 return [] end
    d == 1 && return [-p[0] / p[1]]

    chopped_trimmed = chop(truncate(p))
    n_trail = length(p) - length(chopped_trimmed)
    comp = companion(chopped_trimmed)
    L = eigvals(rot180(comp))
    append!(L, zeros(n_trail))
    return sort!(L, rev = true, by = norm)
end

"""
    companion(::AbstractPolynomial)

Return the companion matrix for the given polynomial.

# References
[Companion Matrix](https://en.wikipedia.org/wiki/Companion_matrix)
"""
companion(::AbstractPolynomial)

"""
    vander(::Type{AbstractPolynomial}, x::AbstractVector, deg::Integer)

Calculate the psuedo-Vandermonde matrix of the given polynomial type with the given degree.

# References
[Vandermonde Matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix)
"""
vander(::Type{<:AbstractPolynomial}, x::AbstractVector, deg::Integer)

"""
    integral(::AbstractPolynomial, k=0)

Returns a polynomial that is the integral of the given polynomial with constant term `k` added.
"""
integral(::AbstractPolynomial, k::Number)
integral(p::AbstractPolynomial) = integral(p, 0)

"""
    integrate(::AbstractPolynomial, a, b)

Compute the definite integral of the given polynomial from `a` to `b`. Will throw an error if either `a` or `b` are out of the polynomial's domain.
"""
function integrate(p::AbstractPolynomial, a::Number, b::Number)
    P = integral(p)
    return P(b) - P(a)
end

"""
    derivative(::AbstractPolynomial, order::Int = 1)

Returns a polynomail that is the `order`th derivative of the given polynomial. `order` must be non-negative.
"""
derivative(::AbstractPolynomial, ::Int)
derivative(p::AbstractPolynomial) = derivative(p, 1)

"""
    truncate!(::AbstractPolynomial{T}; 
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

In-place version of [`truncate`](@ref)
"""
function truncate!(p::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0) where {T}
    max_coeff = maximum(abs, p.coeffs)
    thresh = max_coeff * rtol + atol
    map!(c->abs(c) <= thresh ? zero(T) : c, p.coeffs, p.coeffs)
    return chop!(p, rtol = rtol, atol = atol)
end

"""
    truncate(::AbstractPolynomial{T}; 
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

Rounds off coefficients close to zero, as determined by `rtol` and `atol`, and then chops any leading zeros. Returns a new polynomial.
"""
function Base.truncate(p::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0) where {T}
    truncate!(deepcopy(p), rtol = rtol, atol = atol)
end

"""
    chop!(::AbstractPolynomial{T}; 
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0))

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
    chop(::AbstractPolynomial{T}; 
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0))

Removes any leading coefficients that are approximately 0 (using `rtol` and `atol`). Returns a polynomial whose degree will guaranteed to be equal to or less than the given polynomial's.
"""
function Base.chop(p::AbstractPolynomial{T}; rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0) where {T}
    chop!(deepcopy(p), rtol = rtol, atol = atol)
end

"""
    variable(var=:x)
    variable(::Type{<:AbstractPolynomial}, var=:x)
    variable(p::AbstractPolynomial, var=p.var)

Return the indeterminate of a given polynomial. If no type is give, will default to [`Polynomial`](@ref)

# Examples
```jldoctest
julia> x = variable()
Polynomial(x)

julia> p = 100 + 24x - 3x^2
Polynomial(100 + 24*x - 3*x^2)

julia> roots((x - 3) * (x + 2))
2-element Array{Float64,1}:
  3.0
 -2.0

```
"""
variable(::Type{P}, var::SymbolLike = :x) where P <: AbstractPolynomial = P([0, 1], var)
variable(p::AbstractPolynomial, var::SymbolLike = p.var) = variable(typeof(p), var)
variable(var::SymbolLike = :x) = variable(Polynomial{Int})

#=
Linear Algebra
=#
"""
    norm(::AbstractPolynomial, p=2)

Calculates the p-norm of the polynomial's coefficients
"""
LinearAlgebra.norm(q::AbstractPolynomial, p::Real = 2) = norm(coeffs(q), p)
"""
    conj(::AbstractPolynomial)

Returns the complex conjugate of the polynomial
"""
LinearAlgebra.conj(p::P) where {P <: AbstractPolynomial} = P(conj(coeffs(p)))
LinearAlgebra.transpose(p::AbstractPolynomial) = p
LinearAlgebra.transpose!(p::AbstractPolynomial) = p

#=
Conversions
=#
Base.convert(::Type{P}, p::P) where {P <: AbstractPolynomial} = p
Base.convert(P::Type{<:AbstractPolynomial}, x) = P(x)
Base.promote_rule(::Type{<:AbstractPolynomial{T}}, ::Type{<:AbstractPolynomial{S}}) where {T,S} = Polynomial{promote_type(T, S)}

#=
Inspection
=#
"""
    length(::AbstractPolynomial)

The length of the polynomial.
"""
Base.length(p::AbstractPolynomial) = length(p.coeffs)
"""
    size(::AbstractPolynomial, [i])

Returns the size of the polynomials coefficients, along axis `i` if provided.
"""
Base.size(p::AbstractPolynomial) = size(p.coeffs)
Base.size(p::AbstractPolynomial, i::Integer) = size(p.coeffs, i)
Base.eltype(p::AbstractPolynomial{T}) where {T} = T
Base.eltype(::Type{P}) where {P <: AbstractPolynomial} = P
function Base.iszero(p::AbstractPolynomial)
    if length(p) == 0 return true end
    return length(p) == 1 && p[0] == 0
end

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

"""
    order(::AbstractPolynomial)

The order of the polynomial. This is the same as [`length`](@ref).
"""
order(p::AbstractPolynomial) = length(p)
hasnan(p::AbstractPolynomial) = any(isnan.(p.coeffs))

"""
    domain(::Type{<:AbstractPolynomial})

Returns the domain of the polynomial.
"""
domain(::Type{<:AbstractPolynomial})
domain(::P) where {P <: AbstractPolynomial} = domain(P)

"""
    scale_to_domain(::Type{<:AbstractPolynomial}, x)

Given values of x that are assumed to be unbounded (-∞, ∞), return values rescaled to the domain of the given polynomial.

# Examples
```jldoctest
julia> x = -10:10
-10:10

julia> scale_to_domain(ChebyshevT, x)
-1.0:0.1:1.0

```
"""
function scale_to_domain(P::Type{<:AbstractPolynomial}, x)
    d = domain(P)
    x_scaled = x .* (last(d) - first(d)) / (last(x) - first(x))
end
scale_to_domain(::P, x) where {P <: AbstractPolynomial} = scale_to_domain(P, x)

#=
indexing
=#
Base.firstindex(p::AbstractPolynomial) = 0
Base.lastindex(p::AbstractPolynomial) = length(p) - 1
Base.eachindex(p::AbstractPolynomial) = 0:length(p) - 1
Base.broadcastable(p::AbstractPolynomial) = Ref(p)

# iteration
Base.collect(p::P) where {P <: AbstractPolynomial} = collect(P, p)
Base.iterate(p::AbstractPolynomial) = (p[0] * one(typeof(p)), 1)
function Base.iterate(p::AbstractPolynomial, state)
    state <= length(p) - 1 ? (p[state] * variable(p)^(state), state + 1) : nothing
end

# getindex
function Base.getindex(p::AbstractPolynomial{T}, idx::Int) where {T <: Number}
    idx < 0 && throw(BoundsError(p, idx))
    idx ≥ length(p) && return zero(T)
    return p.coeffs[idx + 1]
end
Base.getindex(p::AbstractPolynomial, idx::Number) = getindex(p, convert(Int, idx))
Base.getindex(p::AbstractPolynomial, indices) = [getindex(p, i) for i in indices]
Base.getindex(p::AbstractPolynomial, ::Colon) = p.coeffs

# setindex
function Base.setindex!(p::AbstractPolynomial, value::Number, idx::Int)
    n = length(p.coeffs)
    if n ≤ idx
        resize!(p.coeffs, idx + 1)
        p.coeffs[n + 1:idx] .= 0
    end
    p.coeffs[idx + 1] = value
    return p
end

Base.setindex!(p::AbstractPolynomial, value::Number, idx::Number) = setindex!(p, value, convert(Int, idx))
Base.setindex!(p::AbstractPolynomial, value::Number, indices) = [setindex!(p, value, i) for i in indices]
Base.setindex!(p::AbstractPolynomial, values, indices) = [setindex!(p, v, i) for (v, i) in zip(values, indices)]
Base.setindex!(p::AbstractPolynomial, value::Number, ::Colon) = setindex!(p, value, eachindex(p))
Base.setindex!(p::AbstractPolynomial, values, ::Colon) = [setindex!(p, v, i) for (v, i) in zip(values, eachindex(p))]

#=
identity
=#
Base.copy(p::P) where {P <: AbstractPolynomial} = P(copy(p.coeffs), p.var)
Base.hash(p::AbstractPolynomial, h::UInt) = hash(p.var, hash(p.coeffs, h))
Base.zero(::Type{P}) where {P <: AbstractPolynomial} = P(zeros(1))
Base.zero(p::P) where {P <: AbstractPolynomial} = zero(P)
Base.one(::Type{P}) where {P <: AbstractPolynomial} = P(ones(1))
Base.one(p::P) where {P <: AbstractPolynomial} = one(P)

#=
arithmetic
=#
Base.:-(p::P) where {P <: AbstractPolynomial} = P(-p.coeffs, p.var)
Base.:+(c::Number, p::AbstractPolynomial) = +(p, c)
Base.:-(p::AbstractPolynomial, c::Number) = +(p, -c)
Base.:-(c::Number, p::AbstractPolynomial, ) = +(-p, c)
Base.:*(c::Number, p::AbstractPolynomial) = *(p, c)

function Base.:*(p::P, c::S) where {P <: AbstractPolynomial,S}
    T = promote_type(P, S)
    return T(p.coeffs .* c, p.var)
end
function Base.:/(p::P, c::S) where {T,P <: AbstractPolynomial{T},S}
    R = promote_type(P, eltype(one(T) / one(S)))
    return R(p.coeffs ./ c, p.var)
end
Base.:-(p1::AbstractPolynomial, p2::AbstractPolynomial) = +(p1, -p2)

function Base.:+(p::P, n::Number) where {P <: AbstractPolynomial}
    p1, p2 = promote(p, n)
    return p1 + p2
end

function Base.:+(p1::P, p2::O) where {P <: AbstractPolynomial,O <: AbstractPolynomial}
    p1, p2 = promote(p1, p2)
    return p1 + p2
end

function Base.:*(p1::P, p2::O) where {P <: AbstractPolynomial,O <: AbstractPolynomial}
    p1, p2 = promote(p1, p2)
    return p1 * p2
end

Base.:^(p::AbstractPolynomial, n::Integer) = Base.power_by_squaring(p, n)

function Base.divrem(num::P, den::O) where {P <: AbstractPolynomial,O <: AbstractPolynomial}
    n, d = promote(num, den)
    return divrem(n, d)
end

"""
    gcd(a::Polynomial, b::Polynomial)

Find the greatest common denominator of two polynomials recursively using
[Euclid's algorithm](http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm).

# Examples

```jldoctest
julia> gcd(fromroots([1, 1, 2]), fromroots([1, 2, 3]))
Polynomial(4.0 - 6.0*x + 2.0*x^2)

```
"""
function Base.gcd(p1::AbstractPolynomial{T}, p2::AbstractPolynomial{S}) where {T, S}
    r₀, r₁ = promote(p1, p2)
    iter    = 1
    itermax = length(r₁)
  
    while r₁ ≉ zero(r₁) && iter ≤ itermax   # just to avoid unnecessary recursion
      _, rtemp  = divrem(r₀, r₁)
      r₀        = r₁
      r₁        = truncate(rtemp)
      iter      += 1
    end
    return r₀
end

"""
    div(::AbstractPolynomial, ::AbstractPolynomial)
"""
Base.div(n::AbstractPolynomial, d::AbstractPolynomial) = divrem(n, d)[1]

"""
    rem(::AbstractPolynomial, ::AbstractPolynomial)
"""
Base.rem(n::AbstractPolynomial, d::AbstractPolynomial) = divrem(n, d)[2]

#=
Comparisons
=#
Base.isequal(p1::P, p2::P) where {P <: AbstractPolynomial} = hash(p1) == hash(p2)
Base.:(==)(p1::AbstractPolynomial, p2::AbstractPolynomial) = (p1.var == p2.var) && (p1.coeffs == p2.coeffs)
Base.:(==)(p::AbstractPolynomial, n::Number) = p.coeffs == [n]
Base.:(==)(n::Number, p::AbstractPolynomial) = p == n

function Base.isapprox(p1::AbstractPolynomial{T}, p2::AbstractPolynomial{S}; rtol::Real = (Base.rtoldefault(T, S, 0)), atol::Real = 0) where {T,S}
    p1, p2 = promote(p1, p2)
    if p1.var != p2.var error("p1 and p2 must have same var") end
    p1t = truncate(p1; rtol = rtol, atol = atol)
    p2t = truncate(p2; rtol = rtol, atol = atol)
    if length(p1t) ≠ length(p2t) return false end
    isapprox(p1t.coeffs, p2t.coeffs, rtol = rtol, atol = atol)
end

function Base.isapprox(p1::AbstractPolynomial{T}, n::S; rtol::Real = (Base.rtoldefault(T, S, 0)), atol::Real = 0) where {T,S}
    p1t = truncate(p1, rtol = rtol, atol = atol)
    if length(p1t) != 1 return false end
    isapprox(p1t.coeffs, [n], rtol = rtol, atol = atol)
end

Base.isapprox(n::S, p1::AbstractPolynomial{T}; rtol::Real = (Base.rtoldefault(T, S, 0)), atol::Real = 0) where {T,S} = isapprox(p1, n, rtol = rtol, atol = atol)
