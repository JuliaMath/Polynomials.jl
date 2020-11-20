using LinearAlgebra

export fromroots,
       truncate!,
       chop!,
       coeffs,
       degree,
       domain,
       mapdomain,
       order,
       hasnan,
       roots,
       companion,
       vander,
       fit,
       integrate,
       derivative,
       variable,
       isintegral,
       ismonic

"""
    fromroots(::AbstractVector{<:Number}; var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractVector{<:Number}; var=:x)

Construct a polynomial of the given type given the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest common
julia> using Polynomials

julia> r = [3, 2]; # (x - 3)(x - 2)

julia> fromroots(r)
Polynomial(6 - 5*x + x^2)
```
"""
function fromroots(P::Type{<:AbstractPolynomial}, roots::AbstractVector; var::SymbolLike = :x)
    x = variable(P, var)
    p =  prod(x - r for r in roots)
    return truncate!(p)
end
fromroots(r::AbstractVector{<:Number}; var::SymbolLike = :x) =
    fromroots(Polynomial, r, var = var)

"""
    fromroots(::AbstractMatrix{<:Number}; var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractMatrix{<:Number}; var=:x)

Construct a polynomial of the given type using the eigenvalues of the given matrix as the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest common
julia> using Polynomials

julia> A = [1 2; 3 4]; # (x - 5.37228)(x + 0.37228)

julia> fromroots(A)
Polynomial(-1.9999999999999998 - 5.0*x + 1.0*x^2)
```
"""
fromroots(P::Type{<:AbstractPolynomial},
    A::AbstractMatrix{T};
    var::SymbolLike = :x,) where {T <: Number} = fromroots(P, eigvals(A), var = var)
fromroots(A::AbstractMatrix{T}; var::SymbolLike = :x) where {T <: Number} =
    fromroots(Polynomial, eigvals(A), var = var)

"""
    fit(x, y, deg=length(x) - 1; [weights], var=:x)
    fit(::Type{<:AbstractPolynomial}, x, y, deg=length(x)-1; [weights], var=:x)

Fit the given data as a polynomial type with the given degree. Uses
linear least squares to minimize the norm of `V⋅c - y`, where `V` is
the Vandermonde matrix and `c` are the coefficients of the polynomial
fit.

This will automatically scale your data to the [`domain`](@ref) of the
polynomial type using [`mapdomain`](@ref). The default polynomial type
is [`Polynomial`](@ref).

When weights are given, as either a `Number`, `Vector` or `Matrix`,
this will use weighted linear least squares. That is, the norm of
`W ⋅ (y - V ⋅ x)` is minimized. (As of now, the weights are specified 
using their squares: for a number use `w^2`, for a vector `wᵢ^2`, and for a matrix
 specify `W'*W`. This behavior may change in the future.)

"""
function fit(P::Type{<:AbstractPolynomial},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::Integer = length(x) - 1;
             weights = nothing,
             var = :x,) where {T}
    _fit(P, x, y, deg; weights=weights, var=var)
end

fit(P::Type{<:AbstractPolynomial},
    x,
    y,
    deg::Integer = length(x) - 1;
    weights = nothing,
    var = :x,) = fit′(P, promote(collect(x), collect(y))..., deg; weights = weights, var = var)

#  avoid issue  214
fit′(P::Type{<:AbstractPolynomial}, x, y, args...;kwargs...) = throw(MethodError("x and y do not produce abstract   vectors"))
fit′(P::Type{<:AbstractPolynomial},
     x::AbstractVector{T},
     y::AbstractVector{T},
     args...; kwargs...) where {T} = fit(P,x,y,args...;  kwargs...)
         
         
fit(x::AbstractVector,
    y::AbstractVector,
    deg::Integer = length(x) - 1;
    weights = nothing,
    var = :x,) = fit(Polynomial, x, y, deg; weights = weights, var = var)

function _fit(P::Type{<:AbstractPolynomial},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::Integer = length(x) - 1;
             weights = nothing,
             var = :x,) where {T}
    x = mapdomain(P, x)
    vand = vander(P, x, deg)
    if weights !== nothing
        coeffs = _wlstsq(vand, y, weights)
    else
        coeffs = qr(vand) \ y
    end
    R = float(T)
    return P(R.(coeffs), var)
end


# Weighted linear least squares
# TODO: Breaking change for 2.0: use non-squared weights
_wlstsq(vand, y, W::Number) = _wlstsq(vand, y, fill!(similar(y), W))
function _wlstsq(vand, y, w::AbstractVector)
    W = Diagonal(sqrt.(w))
    qr(W * vand) \ (W * y)
end
_wlstsq(vand, y, W::AbstractMatrix) = qr(vand' * W * vand) \ (vand' * W * y)

"""
    roots(::AbstractPolynomial; kwargs...)

Returns the roots of the given polynomial. This is calculated via the eigenvalues of the companion matrix. The `kwargs` are passed to the `LinearAlgeebra.eigvals` call.

!!! note

        The [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl) package provides an alternative that is a bit faster and a bit more accurate; the [FastPolynomialRoots](https://github.com/andreasnoack/FastPolynomialRoots.jl) provides an interface to FORTRAN code implementing an algorithm that can handle very large polynomials (it is  `O(n^2)` not `O(n^3)`. The [AMRVW.jl](https://github.com/jverzani/AMRVW.jl) package implements the algorithm in Julia, allowing the use of other  number types.

"""
function roots(q::AbstractPolynomial{T}; kwargs...) where {T <: Number}

    p = convert(Polynomial{T},  q)
    roots(p; kwargs...)

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
    integrate(::AbstractPolynomial, C=0)

Returns the indefinite integral of the polynomial with constant `C`.
"""
integrate(p::AbstractPolynomial, C::Number = 0) = integrate(p, C)

"""
    integrate(::AbstractPolynomial, a, b)

Compute the definite integral of the given polynomial from `a` to `b`. Will throw an error if either `a` or `b` are out of the polynomial's domain.
"""
function integrate(p::AbstractPolynomial, a::Number, b::Number)
    P = integrate(p)
    return P(b) - P(a)
end

"""
    derivative(::AbstractPolynomial, order::Int = 1)

Returns a polynomial that is the `order`th derivative of the given polynomial. `order` must be non-negative.
"""
derivative(::AbstractPolynomial, ::Int)

"""
    truncate!(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

In-place version of [`truncate`](@ref)
"""
function truncate!(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {T}
    max_coeff = maximum(abs, coeffs(p))
    thresh = max_coeff * rtol + atol
    map!(c->abs(c) <= thresh ? zero(T) : c, coeffs(p), coeffs(p))
    return chop!(p, rtol = rtol, atol = atol)
end

"""
    truncate(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

Rounds off coefficients close to zero, as determined by `rtol` and `atol`, and then chops any leading zeros. Returns a new polynomial.
"""
function Base.truncate(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
    atol::Real = 0,) where {T}
    truncate!(deepcopy(p), rtol = rtol, atol = atol)
end

"""
    chop!(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0))

In-place version of [`chop`](@ref)
"""
function chop!(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
               atol::Real = 0,) where {T}
    isempty(values(p)) && return p
    tol = norm(p) * rtol + atol
    for i = lastindex(p):-1:0
        val = p[i]
        if abs(val) > tol #!isapprox(val, zero(T); rtol = rtol, atol = atol)
            resize!(p.coeffs, i + 1); 
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
function Base.chop(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
    atol::Real = 0,) where {T}
    chop!(deepcopy(p), rtol = rtol, atol = atol)
end




"""
    check_same_variable(p::AbstractPolynomial, q::AbstractPolynomial)

Check if either `p` or `q` is constant or if `p` and `q` share the same variable
"""
check_same_variable(p::AbstractPolynomial, q::AbstractPolynomial) =
    (Polynomials.isconstant(p) || Polynomials.isconstant(q)) || p.var ==  q.var

#=
Linear Algebra =#
"""
    norm(::AbstractPolynomial, p=2)

Calculates the p-norm of the polynomial's coefficients
"""
function LinearAlgebra.norm(q::AbstractPolynomial, p::Real = 2)
    vs = values(q)
    return norm(vs, p) # if vs=() must be handled in special type
end

"""
    conj(::AbstractPolynomial)

Returns the complex conjugate of the polynomial
"""
LinearAlgebra.conj(p::P) where {P <: AbstractPolynomial} = map(conj, p)
LinearAlgebra.adjoint(p::P) where {P <: AbstractPolynomial} = map(adjoint, p) 
LinearAlgebra.transpose(p::AbstractPolynomial) = p
LinearAlgebra.transpose!(p::AbstractPolynomial) = p

#=
Conversions =#
Base.convert(::Type{P}, p::P) where {P <: AbstractPolynomial} = p
Base.convert(P::Type{<:AbstractPolynomial}, x) = P(x)
Base.promote_rule(::Type{<:AbstractPolynomial{T}},
    ::Type{<:AbstractPolynomial{S}},
) where {T,S} = Polynomial{promote_type(T, S)}

#=
Inspection =#
"""
    length(::AbstractPolynomial)

The length of the polynomial.
"""
Base.length(p::AbstractPolynomial) = length(coeffs(p))
"""
    size(::AbstractPolynomial, [i])

Returns the size of the polynomials coefficients, along axis `i` if provided.
"""
Base.size(p::AbstractPolynomial) = size(coeffs(p))
Base.size(p::AbstractPolynomial, i::Integer) = size(coeffs(p), i)
Base.eltype(p::AbstractPolynomial{T}) where {T} = T
# in  analogy  with  polynomial as a Vector{T} with different operations defined.
Base.eltype(::Type{<:AbstractPolynomial}) = Float64
Base.eltype(::Type{<:AbstractPolynomial{T}}) where {T} = T
#Base.eltype(::Type{P}) where {P <: AbstractPolynomial} = P # changed  in v1.1.0
Base.iszero(p::AbstractPolynomial) = all(iszero, p)

# See discussions in https://github.com/JuliaMath/Polynomials.jl/issues/258
"""
    all(pred, poly::AbstractPolynomial)

Test whether all coefficients of an `AbstractPolynomial` satisfy predicate `pred`.

You can implement `isreal`, etc., to a `Polynomial` by using `all`.
"""
Base.all(pred, p::AbstractPolynomial) = all(pred, values(p))
"""
    any(pred, poly::AbstractPolynomial)

Test whether any coefficient of an `AbstractPolynomial` satisfies predicate `pred`.
"""
Base.any(pred, p::AbstractPolynomial) = any(pred, values(p))




"""
    map(fn, p::AbstractPolynomial, args...)

Transform coefficients of `p` by applying a function (or other callables) `fn` to each of them.

You can implement `real`, etc., to a `Polynomial` by using `map`.
"""
Base.map(fn, p::P, args...)  where {P<:AbstractPolynomial} = _convert(p, map(fn, coeffs(p), args...))

"""
    isreal(p::AbstractPolynomial)

Determine whether a polynomial is a real polynomial, i.e., having only real numbers as coefficients.

See also: [`real`](@ref)
"""
Base.isreal(p::AbstractPolynomial) = all(isreal, p)
"""
    real(p::AbstractPolynomial)

Construct a real polynomial from the real parts of the coefficients of `p`.

See also: [`isreal`](@ref)

!!! note
    This could cause losing terms in `p`. This method is usually called on polynomials like `p = Polynomial([1, 2 + 0im, 3.0, 4.0 + 0.0im])` where you want to chop the imaginary parts of the coefficients of `p`.
"""
Base.real(p::AbstractPolynomial) = map(real, p)

"""
    isintegral(p::AbstractPolynomial)

Determine whether a polynomial is an integer polynomial, i.e., having only integers as coefficients.
"""
isintegral(p::AbstractPolynomial) = all(isinteger, p)

"""
    ismonic(p::AbstractPolynomial)

Determine whether a polynomial is a monic polynomial, i.e., its leading coefficient is one.
"""
ismonic(p::AbstractPolynomial) = isone(p[end])

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
degree(p::AbstractPolynomial) = iszero(p) ? -1 : lastindex(p) 


"""
    isconstant(::AbstractPolynomial)

Is the polynomial  `p` a constant.
"""
isconstant(p::AbstractPolynomial) = degree(p) <= 0




hasnan(p::AbstractPolynomial) = any(isnan, p)

"""
    domain(::Type{<:AbstractPolynomial})

Returns the domain of the polynomial.
"""
domain(::Type{<:AbstractPolynomial})
domain(::P) where {P <: AbstractPolynomial} = domain(P)

"""
    mapdomain(::Type{<:AbstractPolynomial}, x::AbstractArray)
    mapdomain(::AbstractPolynomial, x::AbstractArray)

Given values of x that are assumed to be unbounded (-∞, ∞), return values rescaled to the domain of the given polynomial.

# Examples
```jldoctest  common
julia> using Polynomials

julia> x = -10:10
-10:10

julia> extrema(mapdomain(ChebyshevT, x))
(-1.0, 1.0)

```
"""
function mapdomain(P::Type{<:AbstractPolynomial}, x::AbstractArray)
    d = domain(P)
    x = collect(x)
    x_zerod = x .- minimum(x)
    x_scaled = x_zerod .* (last(d) - first(d)) ./ maximum(x_zerod)
    x_scaled .+= first(d)
    return x_scaled
end
mapdomain(::P, x::AbstractArray) where {P <: AbstractPolynomial} = mapdomain(P, x)
#=
indexing =#
Base.firstindex(p::AbstractPolynomial) = 0
Base.lastindex(p::AbstractPolynomial) = length(p) - 1
Base.eachindex(p::AbstractPolynomial) = 0:length(p) - 1
Base.broadcastable(p::AbstractPolynomial) = Ref(p)
# like coeffs, though possibly only non-zero values (e.g. SparsePolynomial)
Base.values(p::AbstractPolynomial) = coeffs(p) 

# iteration
# iteration occurs over the basis polynomials
Base.iterate(p::AbstractPolynomial) = (p[0] * one(typeof(p)), 1)
function Base.iterate(p::AbstractPolynomial, state)
    state <= length(p) - 1 ? (p[state] * basis(p, state), state + 1) : nothing
end


Base.collect(p::P) where {P <: AbstractPolynomial} = collect(P, p)

# getindex
function Base.getindex(p::AbstractPolynomial{T}, idx::Int) where {T <: Number}
    idx < 0 && throw(BoundsError(p, idx))
    idx ≥ length(p) && return zero(T)
    return coeffs(p)[idx + 1]
end
Base.getindex(p::AbstractPolynomial, idx::Number) = getindex(p, convert(Int, idx))
Base.getindex(p::AbstractPolynomial, indices) = [getindex(p, i) for i in indices]
Base.getindex(p::AbstractPolynomial, ::Colon) = coeffs(p)

# setindex
function Base.setindex!(p::AbstractPolynomial, value::Number, idx::Int)
    n = length(coeffs(p))
    if n ≤ idx
        resize!(p.coeffs, idx + 1)
        p.coeffs[n + 1:idx] .= 0
    end
    p.coeffs[idx + 1] = value
    return p
end

Base.setindex!(p::AbstractPolynomial, value::Number, idx::Number) =
    setindex!(p, value, convert(Int, idx))
Base.setindex!(p::AbstractPolynomial, value::Number, indices) =
    [setindex!(p, value, i) for i in indices]
Base.setindex!(p::AbstractPolynomial, values, indices) =
    [setindex!(p, v, i) for (v, i) in zip(values, indices)]
Base.setindex!(p::AbstractPolynomial, value::Number, ::Colon) =
    setindex!(p, value, eachindex(p))
Base.setindex!(p::AbstractPolynomial, values, ::Colon) =
    [setindex!(p, v, i) for (v, i) in zip(values, eachindex(p))]

#=
identity =#
Base.copy(p::P) where {P <: AbstractPolynomial} = _convert(p, copy(coeffs(p)))
Base.hash(p::AbstractPolynomial, h::UInt) = hash(p.var, hash(coeffs(p), h))

#=
zero, one, variable, basis =#
"""
    zero(::Type{<:AbstractPolynomial})
    zero(::AbstractPolynomial)

Returns a representation of 0 as the given polynomial.
"""
Base.zero(::Type{P}, var=:x) where {P <: AbstractPolynomial} = ⟒(P)(zeros(eltype(P), 1), var)
Base.zero(p::P) where {P <: AbstractPolynomial} = zero(P, p.var)
"""
    one(::Type{<:AbstractPolynomial})
    one(::AbstractPolynomial)

Returns a representation of 1 as the given polynomial.
"""
Base.one(::Type{P}, var=:x) where {P <: AbstractPolynomial} = ⟒(P)(ones(eltype(P),1), var)  # assumes  p₀ = 1
Base.one(p::P) where {P <: AbstractPolynomial} = one(P, p.var)

Base.oneunit(::Type{P}, args...) where {P <: AbstractPolynomial} = one(P, args...)
Base.oneunit(p::P, args...) where {P <: AbstractPolynomial} = one(p, args...)


"""
    variable(var=:x)
    variable(::Type{<:AbstractPolynomial}, var=:x)
    variable(p::AbstractPolynomial, var=p.var)

Return the monomial `x` in the indicated polynomial basis.  If no type is give, will default to [`Polynomial`](@ref). Equivalent  to  `P(var)`.

# Examples
```jldoctest  common
julia> using Polynomials

julia> x = variable()
Polynomial(x)

julia> p = 100 + 24x - 3x^2
Polynomial(100 + 24*x - 3*x^2)

julia> roots((x - 3) * (x + 2))
2-element Array{Float64,1}:
 -2.0
  3.0

```
"""
variable(::Type{P}, var::SymbolLike = :x) where {P <: AbstractPolynomial} = MethodError()
variable(p::AbstractPolynomial, var::SymbolLike = p.var) = variable(typeof(p), var)
variable(var::SymbolLike = :x) = variable(Polynomial{Int}, var)

# basis
# var is a positional argument, not a keyword; can't deprecate so we do `_var; var=_var`
#@deprecate basis(p::P, k::Int; var=:x)  where {P<:AbstractPolynomial}  basis(p, k, var)
#@deprecate basis(::Type{P}, k::Int; var=:x) where {P <: AbstractPolynomial} basis(P, k,var)
# return the kth basis polynomial for the given polynomial type, e.g. x^k for Polynomial{T}
function basis(::Type{P}, k::Int, _var::SymbolLike=:x; var=_var) where {P <: AbstractPolynomial}
    zs = zeros(Int, k+1)
    zs[end] = 1
    ⟒(P){eltype(P)}(zs, var)
end
basis(p::P, k::Int, _var::SymbolLike=:x; var=_var) where {P<:AbstractPolynomial} = basis(P, k, var)

#=
arithmetic =#
Base.:-(p::P) where {P <: AbstractPolynomial} = _convert(p, -coeffs(p)) 
Base.:+(c::Number, p::AbstractPolynomial) = +(p, c)
Base.:-(p::AbstractPolynomial, c::Number) = +(p, -c)
Base.:-(c::Number, p::AbstractPolynomial) = +(-p, c)
Base.:*(c::Number, p::AbstractPolynomial) = *(p, c)

function Base.:*(p::P, c::S) where {P <: AbstractPolynomial,S}
    _convert(p, coeffs(p) .* c)
end

function Base.:/(p::P, c::S) where {P <: AbstractPolynomial,S}
    _convert(p, coeffs(p) ./ c)
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
    gcd(a::AbstractPolynomial, b::AbstractPolynomial; atol::Real=0, rtol::Real=Base.rtoldefault)

Find the greatest common denominator of two polynomials recursively using
[Euclid's algorithm](http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm).

# Examples

```jldoctest common
julia> using Polynomials

julia> gcd(fromroots([1, 1, 2]), fromroots([1, 2, 3]))
Polynomial(4.0 - 6.0*x + 2.0*x^2)

```
"""
function Base.gcd(p1::AbstractPolynomial{T}, p2::AbstractPolynomial{S}; kwargs...) where {T,S}
    gcd(promote(p1, p2)...; kwargs...)
end

function Base.gcd(p1::AbstractPolynomial{T}, p2::AbstractPolynomial{T};
                  atol::Real=zero(real(T)),
                  rtol::Real=Base.rtoldefault(real(T))
                  ) where {T}


    r₀, r₁ = p1, p2
    iter = 1
    itermax = length(r₁)

    while !iszero(r₁) && iter ≤ itermax
        _, rtemp = divrem(r₀, r₁)
        r₀ = r₁
        r₁ = truncate(rtemp; atol=atol, rtol=rtol)  
        iter += 1
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
Comparisons =#
Base.isequal(p1::P, p2::P) where {P <: AbstractPolynomial} = hash(p1) == hash(p2)
Base.:(==)(p1::AbstractPolynomial, p2::AbstractPolynomial) =
    check_same_variable(p1,p2) && (coeffs(p1) == coeffs(p2))
Base.:(==)(p::AbstractPolynomial, n::Number) = degree(p) <= 0 && p[0] == n
Base.:(==)(n::Number, p::AbstractPolynomial) = p == n

function Base.isapprox(p1::AbstractPolynomial{T},
    p2::AbstractPolynomial{S};
    rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {T,S}
    
    p1, p2 = promote(p1, p2)
    check_same_variable(p1, p2)  || error("p1 and p2 must have same var")
    # copy over from abstractarray.jl
    Δ  = norm(p1-p2)
    if isfinite(Δ)
        return Δ <= max(atol, rtol*max(norm(p1), norm(p2)))
    else
        for i in 0:max(degree(p1), degree(p2))
            isapprox(p1[i], p2[i]; rtol=rtol, atol=atol) || return false
        end
        return true
    end
end

function Base.isapprox(p1::AbstractPolynomial{T},
                       n::S;
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {T,S}
    return isapprox(p1, _convert(p1, [n])) 
end

Base.isapprox(n::S,
    p1::AbstractPolynomial{T};
    rtol::Real = (Base.rtoldefault(T, S, 0)),
    atol::Real = 0,) where {T,S} = isapprox(p1, n, rtol = rtol, atol = atol)
