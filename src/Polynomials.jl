# Poly type manipulations

isdefined(Base, :__precompile__) && __precompile__()

module Polynomials
#todo: sparse polynomials?

using Compat

export Poly, poly
export degree, coeffs, variable
export polyval, polyint, polyder, roots, polyfit
export Pade, padeval

import Base: start, next, done, length, size, eltype
import Base: endof, getindex, setindex!, copy, zero, one, convert, norm, gcd
import Base: show, print, *, /, //, -, +, ==, isapprox, divrem, div, rem, eltype
import Base: promote_rule, truncate, chop,  conj, transpose, dot, hash
import Base: isequal

const SymbolLike = Union{AbstractString,Char,Symbol}

"""
    Poly{T<:Number}(a::AbstractVector{T}, [x])

Construct a polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct this through
`Poly([a_0, a_1, ..., a_n])`.

The usual arithmetic operators are overloaded to work with polynomials as well as
with combinations of polynomials and scalars. However, operations involving two
polynomials of different variables causes an error.

# Examples

```julia
julia> Poly([1, 0, 3, 4])
Poly(1 + 3⋅x^2 + 4⋅x^3)

julia> Poly([1, 2, 3], :s)
Poly(1 + 2⋅s + 3⋅s^2)

julia> a = Poly([1, 2, 3], :x); b = Poly([1, 2, 3], :s);

julia> a + b
ERROR: Polynomials must have same variable
...

julia> p = Poly([1, 2])
Poly(1 + 2⋅x)

julia> q = Poly([1, 0, -1])
Poly(1 - x^2)

julia> 2p
Poly(2 + 4⋅x)

julia> 2 + p
Poly(3 + 2⋅x)

julia> p - q
Poly(2⋅x + x^2)

julia> p * q
Poly(1 + 2⋅x - x^2 - 2⋅x^3)

julia> q / 2
Poly(0.5 - 0.5⋅x^2)
```
"""
immutable Poly{T}
  a::Vector{T}
  var::Symbol
  @compat function (::Type{Poly}){T<:Number}(a::AbstractVector{T}, var::SymbolLike = :x)
    # if a == [] we replace it with a = [0]
    if length(a) == 0
      return new{T}(zeros(T,1), @compat Symbol(var))
    else
      # determine the last nonzero element and truncate a accordingly
      a_last = max(1,findlast(x->x!=zero(T), a))
      new{T}(a[1:a_last], @compat Symbol(var))
    end
  end
end

Poly(n::Number, var::SymbolLike = :x) = Poly([n], var)
@compat (::Type{Poly{T}}){T,S}(x::AbstractVector{S}, var::SymbolLike = :x) =
  Poly(convert(Vector{T}, x), var)

# create a Poly object from its roots
"""
    poly(r)

Construct a polynomial from its roots. Compare this to the `Poly` type
constructor, which constructs a polynomial from its coefficients.

If `r` is a vector, the constructed polynomial is
``(x - r_1) (x - r_2) \\cdots (x - r_n)``.
If `r` is a matrix, the constructed polynomial is
``(x - e_1) \\cdots (x - e_n)``, where ``e_i`` is the ``i``th eigenvalue
of `r`.

# Examples

```julia
julia> poly([1, 2, 3])   # The polynomial (x - 1)(x - 2)(x - 3)
Poly(-6 + 11⋅x - 6⋅x^2 + x^3)

julia> poly([1 2; 3 4])  # The polynomial (x - 5.37228)(x + 0.37228)
Poly(-1.9999999999999998 - 5.0⋅x + 1.0⋅x^2)
```
"""
function poly{T}(r::AbstractVector{T}, var::SymbolLike=:x)
    n = length(r)
    c = zeros(T, n+1)
    c[1] = one(T)
    for j = 1:n
        for i = j:-1:1
            c[i+1] = c[i+1]-r[j]*c[i]
        end
    end
    return Poly(reverse(c), var)
end
poly(A::Matrix, var::SymbolLike=:x) = poly(eigvals(A), var)


include("show.jl") # display polynomials.

convert{T}(::Type{Poly{T}}, p::Poly{T}) = p
convert{T}(::Type{Poly{T}}, p::Poly) = Poly(convert(Vector{T}, p.a), p.var)
convert{T, S<:Number}(::Type{Poly{T}}, x::S, var::SymbolLike=:x) = Poly(T[x], var)
convert{T, S<:Number}(::Type{Poly{T}}, x::AbstractArray{S}, var::SymbolLike=:x) = map(el->Poly(T[el],var), x)
promote_rule{T, S}(::Type{Poly{T}}, ::Type{Poly{S}}) = Poly{promote_type(T, S)}
promote_rule{T, S<:Number}(::Type{Poly{T}}, ::Type{S}) = Poly{promote_type(T, S)}
eltype{T}(::Poly{T}) = Poly{T}

length(p::Poly) = length(coeffs(p))
endof(p::Poly)  = length(p) - 1

start(p::Poly)        = start(coeffs(p)) - 1
next(p::Poly, state)  = (temp = zeros(coeffs(p)); temp[state+1] = p[state]; (Poly(temp), state+1))
done(p::Poly, state)  = state > degree(p)
eltype{T}(::Type{Poly{T}}) = Poly{T}

size(p::Poly) = size(p.a)
size(p::Poly, i::Integer) = size(p.a, i)

"""
    degree(p::Poly)

Return the degree of the polynomial `p`, i.e. the highest exponent in the polynomial that
has a nonzero coefficient.
"""
degree(p::Poly) = length(p) - 1

"""
    coeffs(p::Poly)

Return the coefficient vector `[a_0, a_1, ..., a_n]` of a polynomial `p`.
"""
coeffs(p::Poly) = p.a

"""
    variable(p::Poly)
    variable([T::Type,] var)
    variable()

Return the indeterminate of a polynomial, i.e. its variable, as a `Poly` object.
When passed no arguments, this is equivalent to `variable(Float64, :x)`.

# Examples

```julia
julia> variable(Poly([1, 2], :x))
Poly(x)

julia> variable(:y)
Poly(1.0⋅y)

julia> variable()
Poly(1.0⋅x)

julia> variable(Float32, :x)
Poly(1.0f0⋅x)
```
"""
variable{T<:Number}(::Type{T}, var::SymbolLike=:x) = Poly([zero(T), one(T)], var)
variable{T}(p::Poly{T}) = variable(T, p.var)
variable(var::SymbolLike=:x) = variable(Float64, var)

"""
    truncate{T}(p::Poly{T}; reltol::Real = Base.rtoldefault(real(T)), abstol::Real = 0)

Return a polynomial with coefficients `a_i` truncated to zero if `|a_i| <= reltol*maxabs(a)+abstol`.
"""
function truncate{T}(p::Poly{T}; reltol::Real = Base.rtoldefault(real(T)),
  abstol::Real = 0)
    a = coeffs(p)
    amax = maximum(abs,a)
    thresh = amax * reltol + abstol
    anew = map(ai -> abs(ai) <= thresh ? zero(T) : ai, a)
    return Poly(anew, p.var)
end

"""
    chop(p::Poly{T}; reltol::Real = Base.rtoldefault(real(T)), abstol::Real = 0)

Chop off leading values from a polynomial which are approximately zero. The tolerances
`reltol` and `abstol` are passed to `isapprox` to check for zeros.
"""
function chop{T}(p::Poly{T}; reltol::Real = Base.rtoldefault(real(T)),
  abstol::Real = 0)
    c = copy(p.a)
    for k=length(c):-1:1
        if !isapprox(c[k], zero(T); rtol=reltol, atol=abstol)
            resize!(c, k)
            return Poly(c, p.var)
        end
    end

    resize!(c,0)
    Poly(c, p.var)
end

"""
    norm(q::Poly, [p])

Return the `p`-norm of a polynomial `q`. If ``q = q_0 + \\ldots + q_n x^n``, then
the `p`-norm is

``
||q||_p = (|q_0|^p + \\ldots + |q_n|^p)^{1/p}
``
"""
norm(q::Poly, args...) = norm(coeffs(q), args...)


"""
    conj(p::Poly)

Conjugate each coefficient of `p`.
"""
conj{T<:Complex}(p::Poly{T}) = Poly(conj(coeffs(p)))

# Define the no-op `transpose` explicitly to avoid future warnings in Julia
transpose(p::Poly) = p

"""
    getindex(p::Poly, i)

If ``p = a_n x^n + a_{n-1}x^{n-1} + \\ldots + a_1 x^1 + a_0``, then `p[i]` returns ``a_i``.
"""
getindex{T}(p::Poly{T}, idx::Int)                     = (idx ≥ length(p.a) ? zero(T) : p.a[idx+1])
getindex{T}(p::Poly{T}, indices::AbstractVector{Int}) = map(idx->p[idx], indices)
getindex{T}(p::Poly{T}, ::Colon)                      = p[0:length(p)-1]

function setindex!(p::Poly, value, idx::Int)
    n = length(p.a)
    if n ≤ idx
        resize!(p.a, idx+1)
        p.a[n+1:idx] = 0
    end
    p.a[idx+1] = value
    return p
end

function setindex!(p::Poly, values::AbstractVector, indices::AbstractVector{Int})
    for (idx, value) in zip(indices, values)
      setindex!(p, value, idx)
    end
    return p
end

function setindex!(p::Poly, value, indices::AbstractVector{Int})
  for idx in indices
    setindex!(p, value, idx)
  end
  return p
end

setindex!(p::Poly, values, ::Colon) = setindex!(p, values, 0:length(p)-1)

eachindex{T}(p::Poly{T}) = 0:(length(p)-1)

copy(p::Poly) = Poly(copy(p.a), p.var)

zero{T}(p::Poly{T}) = Poly(T[], p.var)
zero{T}(::Type{Poly{T}}) = Poly(T[])
one{T}(p::Poly{T}) = Poly([one(T)], p.var)
one{T}(::Type{Poly{T}}) = Poly([one(T)])

## Overload arithmetic operators for polynomial operations between polynomials and scalars
*{T<:Number,S}(c::T, p::Poly{S}) = Poly(c * p.a, p.var)
*{T<:Number,S}(p::Poly{S}, c::T) = Poly(p.a * c, p.var)
/(p::Poly, c::Number) = Poly(p.a / c, p.var)
-(p::Poly) = Poly(-p.a, p.var)
-{T<:Number}(p::Poly, c::T) = +(p, -c)
+{T<:Number}(c::T, p::Poly) = +(p, c)
function +{S,T<:Number}(p::Poly{S}, c::T)
    U = promote_type(S,T)
    p2 = U == S ? copy(p) : convert(Poly{U}, p)
    p2[0] += c
    return p2
end
function -{T<:Number,S}(c::T, p::Poly{S})
    U = promote_type(S,T)
    p2 = convert(Poly{U}, -p)
    p2[0] += c
    return p2
end

function +{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    Poly([p1[i] + p2[i] for i = 0:max(length(p1),length(p2))], p1.var)
end
function -{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    Poly([p1[i] - p2[i] for i = 0:max(length(p1),length(p2))], p1.var)
end


function *{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)-1
    m = length(p2)-1
    a = zeros(R,m+n+1)

    for i = 0:n
        for j = 0:m
            a[i+j+1] += p1[i] * p2[j]
        end
    end
    Poly(a,p1.var)
end

## older . operators, hack to avoid warning on v0.6
dot_operators = quote
    @compat Base.:.+{T<:Number}(c::T, p::Poly) = +(p, c)
    @compat Base.:.+{T<:Number}(p::Poly, c::T) = +(p, c)
    @compat Base.:.-{T<:Number}(p::Poly, c::T) = +(p, -c)
    @compat Base.:.-{T<:Number}(c::T, p::Poly) = +(p, -c)
    @compat Base.:.*{T<:Number,S}(c::T, p::Poly{S}) = Poly(c * p.a, p.var)
    @compat Base.:.*{T<:Number,S}(p::Poly{S}, c::T) = Poly(p.a * c, p.var)
end
VERSION < v"0.6.0-dev" && eval(dot_operators)


# are any values NaN
hasnan(p::Poly) = reduce(|, (@compat isnan.(p.a)))

function divrem{T, S}(num::Poly{T}, den::Poly{S})
    if num.var != den.var
        error("Polynomials must have same variable")
    end
    m = length(den)-1
    if m == 0 && den[0] == 0
        throw(DivideError())
    end
    R = typeof(one(T)/one(S))
    n = length(num)-1
    deg = n-m+1
    if deg <= 0
        return convert(Poly{R}, zero(num)), convert(Poly{R}, num)
    end

    aQ = zeros(R, deg)
    # aR = deepcopy(num.a)
    # @show num.a
    aR = R[ num.a[i] for i = 1:n+1 ]
    for i = n:-1:m
        quot = aR[i+1] / den[m]
        aQ[i-m+1] = quot
        for j = 0:m
            elem = den[j]*quot
            aR[i-(m-j)+1] -= elem
        end
    end
    pQ = Poly(aQ, num.var)
    pR = Poly(aR, num.var)

    return pQ, pR
end

div(num::Poly, den::Poly) = divrem(num, den)[1]
rem(num::Poly, den::Poly) = divrem(num, den)[2]

==(p1::Poly, p2::Poly) = (p1.var == p2.var && p1.a == p2.a)
==(p1::Poly, n::Number) = (coeffs(p1) == [n])
==(n::Number, p1::Poly) = (p1 == n)

"""
    isapprox{T,S}(p1::Poly{T}, p2::Poly{S}; reltol::Real = Base.rtoldefault(T,S), abstol::Real = 0, norm::Function = vecnorm)

Truncate polynomials `p1` and `p2`, and compare the coefficient vectors using the
given `norm` function. The tolerances `reltol` and `abstol` are passed to both
`truncate` and `isapprox`.
"""
function isapprox{T,S}(p1::Poly{T}, p2::Poly{S};
  reltol::Real = Base.rtoldefault(T,S), abstol::Real = 0, norm::Function = vecnorm)
  p1.var == p2.var || error("Polynomials must have same variable")
  p1t = truncate(p1; reltol = reltol, abstol = abstol)
  p2t = truncate(p2; reltol = reltol, abstol = abstol)
  length(p1t) == length(p2t) && isapprox(coeffs(p1t), coeffs(p2t); rtol = reltol,
    atol = abstol, norm = norm)
end

function isapprox{T,S<:Number}(p1::Poly{T}, n::S; reltol::Real = Base.rtoldefault(T,S),
  abstol::Real = 0)
  p1t = truncate(p1; reltol = reltol, abstol = abstol)
  degree(p1t) == 0 && isapprox(coeffs(p1), [n]; rtol = reltol, atol = abstol)
end

isapprox{T,S<:Number}(n::S, p1::Poly{T}; reltol::Real = Base.rtoldefault(T,S),
  abstol::Real = 0) = isapprox(p1, n; reltol = reltol, abstol = abstol)

hash(f::Poly, h::UInt) = hash(f.var, hash(f.a, h))
isequal(p1::Poly, p2::Poly) = hash(p1) == hash(p2)

"""
    polyval(p::Poly, x::Number)

Evaluate the polynomial `p` at `x` using Horner's method. `Poly` objects
are callable, using this function.

# Examples

```julia
julia> p = Poly([1, 0, -1])
Poly(1 - x^2)

julia> polyval(p, 1)
0

julia> p(1)
0
```
"""
function polyval{T,S}(p::Poly{T}, x::S)
    R = promote_type(T,S)

    lenp = length(p)
    if lenp == 0
        return zero(R) * x
    else
        y = convert(R, p[end])
        for i = (endof(p)-1):-1:0
            y = p[i] + x*y
        end
        return y
    end
end

polyval(p::Poly, v::AbstractArray) = map(x->polyval(p, x), v)

@compat (p::Poly)(x) = polyval(p, x)

"""
    polyint(p::Poly, k::Number=0)

Integrate the polynomial `p` term by term, optionally adding a constant term
`k`. The order of the resulting polynomial is one higher than the order of `p`.

# Examples

```julia
julia> polyint(Poly([1, 0, -1]))
Poly(1.0⋅x - 0.3333333333333333⋅x^3)

julia> polyint(Poly([1, 0, -1]), 2)
Poly(2.0 + 1.0⋅x - 0.3333333333333333⋅x^3)
```
"""
# if we do not have any initial condition, assume k = zero(Int)
polyint{T}(p::Poly{T}) = polyint(p, 0)

# if we have coefficients that have `NaN` representation
function polyint{T<:Union{Real,Complex},S<:Number}(p::Poly{T}, k::S)
    hasnan(p) && return Poly(promote_type(T,S)[NaN])
    _polyint(p, k)
end

# if we have initial condition that can represent `NaN`
function polyint{T,S<:Union{Real,Complex}}(p::Poly{T}, k::S)
    isnan(k) && return Poly(promote_type(T,S)[NaN])
    _polyint(p, k)
end

# if we have both coefficients and initial condition that can take `NaN`
function polyint{T<:Union{Real,Complex},S<:Union{Real,Complex}}(p::Poly{T}, k::S)
    (hasnan(p) || isnan(k)) && return Poly(promote_type(T,S)[NaN])
    _polyint(p, k)
end

# otherwise, catch all
polyint{T,S<:Number}(p::Poly{T}, k::S) = _polyint(p, k)

function _polyint{T,S<:Number}(p::Poly{T}, k::S)
    n = length(p)
    R = promote_type(typeof(one(T)/1), S)
    a2 = Vector{R}(n+1)
    a2[1] = k
    for i = 1:n
        a2[i+1] = p[i-1] / i
    end
    return Poly(a2, p.var)
end

"""
    polyint(p::Poly, a::Number, b::Number)

Compute the definite integral of the polynomial `p` over the interval `[a,b]`.

# Examples

```julia
julia> polyint(Poly([1, 0, -1]), 0, 1)
0.6666666666666667
```
"""
function polyint(p::Poly, a::Number, b::Number)
    P = polyint(p)
    P(b) - P(a)
end

"""
    polyder(p::Poly, k=1)

Compute the `k`th derivative of the polynomial `p`.

# Examples

```julia
julia> polyder(Poly([1, 3, -1]))
Poly(3 - 2⋅x)

julia> polyder(Poly([1, 3, -1]), 2)
Poly(-2)
```
"""
# if we have coefficients that can represent `NaN`s
function polyder{T<:Union{Real,Complex}}(p::Poly{T}, order::Int=1)
    n = length(p)
    order < 0       && error("Order of derivative must be non-negative")
    order == 0      && return p
    hasnan(p)       && return Poly(T[NaN], p.var)
    n <= order      && return Poly(T[], p.var)
    _polyder(p, order)
end

# otherwise
function polyder{T}(p::Poly{T}, order::Int=1)
  n = length(p)
  order < 0   && error("Order of derivative must be non-negative")
  order == 0  && return p
  n <= order  && return Poly(T[], p.var)
  _polyder(p, order)
end

function _polyder{T}(p::Poly{T}, order::Int=1)
  n = length(p)
  a2 = Vector{T}(n-order)
  for i = order:n-1
    a2[i-order+1] = p[i] * prod((i-order+1):i)
  end

  return Poly(a2, p.var)
end

polyint{T}(a::AbstractArray{Poly{T}}, k::Number  = 0) = map(p->polyint(p,k),    a)
polyder{T}(a::AbstractArray{Poly{T}}, order::Int = 1) = map(p->polyder(p,order),a)

##################################################
##
## Some functions on polynomials...


# compute the roots of a polynomial
"""
    roots(p::Poly)

Return the roots (zeros) of `p`, with multiplicity. The number of roots
returned is equal to the order of `p`. The returned roots may be real or
complex.

# Examples

```julia
julia> roots(Poly([1, 0, -1]))
2-element Array{Float64,1}:
 -1.0
  1.0

julia> roots(Poly([1, 0, 1]))
2-element Array{Complex{Float64},1}:
 0.0+1.0im
 0.0-1.0im

julia> roots(Poly([0, 0, 1]))
2-element Array{Float64,1}:
 0.0
 0.0

julia> roots(poly([1,2,3,4]))
4-element Array{Float64,1}:
 4.0
 3.0
 2.0
 1.0
```
"""
function roots{T}(p::Poly{T})
    R = promote_type(T, Float64)
    length(p) == 0 && return zeros(R, 0)
    p = truncate(p)
    num_leading_zeros = 0
    while p[num_leading_zeros] ≈ zero(T)
        if num_leading_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_leading_zeros += 1
    end
    num_trailing_zeros = 0
    while p[end - num_trailing_zeros] ≈ zero(T)
        num_trailing_zeros += 1
    end
    n = endof(p)-(num_leading_zeros + num_trailing_zeros)
    n < 1 && return zeros(R, length(p) - num_trailing_zeros - 1)

    companion = diagm(ones(R, n-1), -1)
    an = p[end-num_trailing_zeros]
    companion[1,:] = -p[(end-num_trailing_zeros-1):-1:num_leading_zeros] / an

    D = eigvals(companion)
    r = zeros(eltype(D),length(p)-num_trailing_zeros-1)
    r[1:n] = D
    return r
end
roots{T}(p::Poly{Rational{T}}) = roots(convert(Poly{promote_type(T, Float64)}, p))

## compute gcd of two polynomials
"""
    gcd(a::Poly, b::Poly)

Find the greatest common denominator of two polynomials recursively using
[Euclid's algorithm](http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm).

# Examples

```julia
julia> gcd(poly([1,1,2]), poly([1,2,3])) # returns (x-1)*(x-2)
Poly(4.0 - 6.0⋅x + 2.0⋅x^2)
```
"""
function gcd{T, S}(a::Poly{T}, b::Poly{S})
  R = typeof(one(T)/one(S))
  degree(b) == 0 && b ≈ zero(b) && return convert(Poly{R}, a)

  s, r = divrem(a, b)
  return gcd(b, r)
end


## Fit degree n polynomial to points
"""
    polyfit(x, y, n=length(x)-1, sym=:x)

Fit a polynomial of degree `n` through the points specified by `x` and `y`,
where `n <= length(x) - 1`, using least squares fit. When `n=length(x)-1`
(the default), the interpolating polynomial is returned. The optional fourth
argument can be used to specify the symbol for the returned polynomial.

# Examples

```julia
julia> xs = linspace(0, pi, 5);

julia> ys = map(sin, xs);

julia> polyfit(xs, ys, 2)
Poly(-0.004902082150108854 + 1.242031920509868⋅x - 0.39535103925413095⋅x^2)
```
"""
function polyfit(x, y, n::Int=length(x)-1, sym::Symbol=:x)
    length(x) == length(y) || throw(DomainError)
    1 <= n <= length(x) - 1 || throw(DomainError)

    #
    # here unsure, whether similar(float(x[1]),...), or similar(x,...)
    # however similar may yield unwanted surprise in case of e.g. x::Int
    #
    A=similar(float(x[1:1]), length(x), n+1)
    #
    # TODO: add support for poly coef bitmap
    # (i.e. polynomial with some terms fixed to zero)
    #
    A[:,1]=1
    for i=1:n
        A[:,i+1]=A[:,i] .* x   # cumulative product more precise than x.^n
    end
    Aqr=qrfact(A)   # returns QR object, not a matrix
    p=Aqr\y         # least squares solution via QR
    Poly(p, sym)
end
polyfit(x,y,sym::Symbol) = polyfit(x,y,length(x)-1, sym)

### Pull in others
include("pade.jl")

end # module Poly
