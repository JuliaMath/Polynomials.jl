# Poly type manipulations

isdefined(Base, :__precompile__) && __precompile__()

module Polynomials
#todo: sparse polynomials?

using Compat

export Poly, poly
export degree, coeffs, variable
export polyval, polyint, polyder, roots, polyfit
export Pade, padeval

import Base: length, endof, getindex, setindex!, copy, zero, one, convert, norm, gcd
import Base: show, print, *, /, //, -, +, ==, divrem, div, rem, eltype, .*, .-, .+
import Base: promote_rule, truncate, chop
if VERSION >= v"0.4"
    import Base.call
end

eps{T}(::Type{T}) = zero(T)
eps{F<:AbstractFloat}(x::Type{F}) = Base.eps(F)
eps{T}(x::Type{Complex{T}}) = eps(T)


"""

* `Poly{T<:Number}(a::Vector)`: Construct a polynomial from its coefficients, lowest order first. That is if `p=a_n x^n + ... + a_2 x^2 + a_1 x^1 + a_0`, we construct this through `Poly([a_0, a_1, ..., a_n])`.

Example:
```
Poly([1,0,3,4])    # Poly(1 + 3x^2 + 4x^3)
```

An optional variable parameter can be added:

```
Poly([1,2,3], :s)       # Poly(1 + 2s + 3s^2)
```

The usual arithmetic operators are overloaded to work on polynomials, and combinations of polynomials and scalars.

```
p = Poly([1,2])        # Poly(1 + 2x)
q = Poly([1, 0, -1])   # Poly(1 - x^2)
2p                     # Poly(2 + 4x)
2+p                    # Poly(3 + 2x)
p - q                  # Poly(2x + x^2)
p*q                    # Poly(1 + 2x - x^2 - 2x^3)
q/2                    # Poly(0.5 - 0.5x^2)
```

Note that operations involving polynomials with different variables will error:

```j
p = Poly([1, 2, 3], :x)
q = Poly([1, 2, 3], :s)
p + q                  # ERROR: Polynomials must have same variable.
```

"""
immutable Poly{T<:Number}
    a::Vector{T}
    var::Symbol
    function Poly(a::Vector{T}, var)
        # if a == [] we replace it with a = [0]
        if length(a) == 0
            return new(zeros(T,1), @compat Symbol(var))
        else
            # determine the last nonzero element and truncate a accordingly
            a_last = max(1,findlast(a))
            new(a[1:a_last], @compat Symbol(var))
        end
    end
end

@compat Poly{T<:Number}(a::Vector{T}, var::Union{AbstractString,Symbol,Char}=:x) = Poly{T}(a, var)

# create a Poly object from its roots
"""

* `poly(r::AbstractVector)`: Construct a polynomial from its
  roots. This is in contrast to the `Poly` constructor, which
  constructs a polynomial from its coefficients.

Example:
```
## Represents (x-1)*(x-2)*(x-3)
poly([1,2,3])     # Poly(-6 + 11x - 6x^2 + x^3)
```
"""
function poly{T}(r::AbstractVector{T}, var=:x)
    n = length(r)
    c = zeros(T, n+1)
    c[1] = 1
    for j = 1:n
        for i = j:-1:1
            c[i+1] = c[i+1]-r[j]*c[i]
        end
    end
    return Poly(reverse(c), var)
end
poly(A::Matrix, var=:x) = poly(eigvals(A), var)
poly(A::Matrix, var::AbstractString) = poly(eigvals(A), @compat Symbol(var))
poly(A::Matrix, var::Char) = poly(eig(A)[1], @compat Symbol(var))


include("show.jl") # display polynomials.

convert{T}(::Type{Poly{T}}, p::Poly{T}) = p
convert{T}(::Type{Poly{T}}, p::Poly) = Poly(convert(Vector{T}, p.a), p.var)
convert{T, S<:Number}(::Type{Poly{T}}, x::S) = Poly(promote_type(T, S)[x])
convert{T, S<:Number,n}(::Type{Poly{T}}, x::Array{S,n}) = map(el->convert(Poly{promote_type(T,S)},el),x)
promote_rule{T, S}(::Type{Poly{T}}, ::Type{Poly{S}}) = Poly{promote_type(T, S)}
promote_rule{T, S<:Number}(::Type{Poly{T}}, ::Type{S}) = Poly{promote_type(T, S)}
eltype{T}(::Poly{T}) = T
eltype{T}(::Type{Poly{T}}) = T

"""

`legnth(p::Poly)`: return length of coefficient vector

"""
length(p::Poly) = length(p.a)
endof(p::Poly) = length(p) - 1

"""

`degree(p::Poly)`: return degree of polynomial `p`

"""
degree(p::Poly) = length(p) - 1

"""

`coeffs(p::Poly)`: return coefficient vector [a_0, a_1, ..., a_n]

"""
coeffs(p::Poly) = p.a

"""

Return the indeterminate of a polynomial, `x`.

* `variable(p::Poly)`: return variable of `p` as a `Poly` object.
* `variable(T<:Number, [:x])`: return poly one(T)*x
* `variable([var::Symbol])`: return polynomial 1x over `Float64`.

"""
variable{T}(p::Poly{T}) = poly(zeros(T,1), p.var)
variable{T<:Number}(::Type{T}, var=:x) = poly(zeros(T,1), var)
variable(var::Symbol=:x) = poly([0.0], var)

"""

`truncate{T}(p::Poly{T}; reltol = eps(T), abstol = eps(T))`: returns a polynomial with coefficients a_i truncated to zero if |a_i| <= reltol*maxabs(a)+abstol

"""
function truncate{T}(p::Poly{T}; reltol = eps(T), abstol = eps(T))
    a = coeffs(p)
    amax = maxabs(a)
    anew = map(ai -> abs(ai) <= amax*reltol+abstol ? zero(T) : ai, a)
    return Poly(anew, p.var)
end

"""

`chop(p::Poly{T}; kwargs...)` chop off leading values which are
approximately zero. The tolerances are passed to `isapprox`.

"""
function chop{T}(p::Poly{T}; reltol=zero(T), abstol=2 * eps(T))
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

* `norm(q::Poly, [p])`: return `p` norm of polynomial `q`

"""
norm(q::Poly, args...) = norm(coeffs(q), args...)


"""

* `conj(p::Poly`): return conjugate of polynomial `p`. (Polynomial with conjugate of each coefficient.)

"""
Base.conj{T<:Complex}(p::Poly{T}) = Poly(conj(coeffs(p)))

"""

* `getindex(p::Poly, i)`: If `p=a_n x^n + a_{n-1}x^{n-1} + ... + a_1 x^1 + a_0`, then `p[i]` returns `a_i`.

"""
getindex{T}(p::Poly{T}, i) = (i+1 > length(p.a) ? zero(T) : p.a[i+1])
getindex{T}(p::Poly{T}, idx::AbstractArray) = map(i->p[i], idx)
function setindex!(p::Poly, v, i)
    n = length(p.a)
    if n < i+1
        resize!(p.a,i+1)
        p.a[n+1:i] = 0
    end
    p.a[i+1] = v
    v
end
function setindex!(p::Poly, vs, idx::AbstractArray)
    [setindex!(p, v, i) for (i,v) in zip(idx, vs)]
    p
end

copy(p::Poly) = Poly(copy(p.a), p.var)

zero{T}(p::Poly{T}) = Poly([zero(T)], p.var)
zero{T}(::Type{Poly{T}}) = Poly(T[])
one{T}(p::Poly{T}) = Poly([one(T)], p.var)
one{T}(::Type{Poly{T}}) = Poly([one(T)])

## Overload arithmetic operators for polynomial operations between polynomials and scalars
*{T<:Number,S}(c::T, p::Poly{S}) = Poly(c * p.a, p.var)
*{T<:Number,S}(p::Poly{S}, c::T) = Poly(p.a * c, p.var)
.*{T<:Number,S}(c::T, p::Poly{S}) = Poly(c * p.a, p.var)
.*{T<:Number,S}(p::Poly{S}, c::T) = Poly(p.a * c, p.var)
/(p::Poly, c::Number) = Poly(p.a / c, p.var)
-(p::Poly) = Poly(-p.a, p.var)
-{T<:Number}(p::Poly, c::T) = +(p, -c)
.-{T<:Number}(p::Poly, c::T) = +(p, -c)
.-{T<:Number}(c::T, p::Poly) = +(p, -c)
+{T<:Number}(c::T, p::Poly) = +(p, c)
.+{T<:Number}(c::T, p::Poly) = +(p, c)
.+{T<:Number}(p::Poly, c::T) = +(p, c)
function +{T<:Number}(p::Poly, c::T)
    if length(p) < 1
        return Poly([c,], p.var)
    else
        p2 = copy(p)
        p2.a[1] += c;
        return p2;
    end
end
function -{T<:Number}(c::T, p::Poly)
    if length(p) < 1
        return Poly(T[c,], p.var)
    else
        p2 = -p;
        p2.a[1] += c;
        return p2;
    end
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

@deprecate /(num::Poly, den::Poly)                div(num::Poly, den::Poly)
div(num::Poly, den::Poly) = divrem(num, den)[1]
rem(num::Poly, den::Poly) = divrem(num, den)[2]

function ==(p1::Poly, p2::Poly)
    if p1.var != p2.var
        return false
    else
        return p1.a == p2.a
    end
end

"""
* `polyval(p::Poly, x::Number)`: Evaluate the polynomial `p` at `x` using Horner's method.

Example:
```
polyval(Poly([1, 0, -1]), 0.1)  # 0.99
```

For `julia` version `0.4` or greater, the `call` method can be used:

```
p = Poly([1,2,3])
p(4)   # 57 = 1 + 2*4 + 3*4^2
```

"""
function polyval{T,S}(p::Poly{T}, x::S)
    R = promote_type(T,S)
    lenp = length(p)
    if lenp == 0
        return zero(R) * x
    else
        y = convert(R, p[end]) + 0*x
        for i = (endof(p)-1):-1:0
            y = p[i] + x*y
        end
        return y
    end
end

polyval(p::Poly, v::AbstractVector) = map(x->polyval(p, x), v)

if VERSION >= v"0.4"
    @compat (p::Poly)(x) = polyval(p, x)
end

"""

* `polyint(p::Poly, k::Number=0)`: Integrate the polynomial `p` term
  by term, optionally adding constant term `k`. The order of the
  resulting polynomial is one higher than the order of `p`.

Examples:
```
polyint(Poly([1, 0, -1]))     # Poly(x - 0.3333333333333333x^3)
polyint(Poly([1, 0, -1]), 2)  # Poly(2.0 + x - 0.3333333333333333x^3)
```

"""
function polyint{T}(p::Poly{T}, k::Number=0)
    n = length(p)
    R = typeof(one(T)/1)
    a2 = Array(R, n+1)
    a2[1] = k
    for i = 1:n
        a2[i+1] = p[i-1] / i
    end
    return Poly(a2, p.var)
end

"""

* `polyder(p::Poly)`: Differentiate the polynomial `p` term by
  term. The order of the resulting polynomial is one lower than the
  order of `p`.

Example:
```
polyder(Poly([1, 3, -1]))   # Poly(3 - 2x)
```
"""
function polyder{T}(p::Poly{T}, order::Int=1)
    n = length(p)
    if order < 0
        error("Order of derivative must be non-negative")
    elseif n <= order
        return Poly(zeros(T,0),p.var)
    elseif order == 0
        return p
    else
        a2 = Array(T, n-order)
        for i = order:n-1
            a2[i-order+1] = p[i] * prod((i-order+1):i)
        end
        return Poly(a2, p.var)
    end
end
Base.ctranspose{T}(p::Poly{T}) = polyder(p)

polyint{T}(a::Array{Poly{T},1}, k::Number  = 0) = [ polyint(p,k) for p in a ]
polyder{T}(a::Array{Poly{T},1}, order::Int = 1) = [ polyder(p,order) for p in a ]
polyint{n,T}(a::Array{Poly{T},n}, k::Number  = 0) = map(p->polyint(p,k),a)
polyder{n,T}(a::Array{Poly{T},n}, order::Int = 1) = map(p->polyder(p,order),a)



##################################################
##
## Some functions on polynomials...


# compute the roots of a polynomial
"""

* `roots(p::Poly)`: Return the roots (zeros) of `p`, with
  multiplicity. The number of roots returned is equal to the order of
  `p`. The returned roots may be real or complex.

Examples:
```
roots(Poly([1, 0, -1]))    # [-1.0, 1.0]
roots(Poly([1, 0, 1]))     # [0.0+1.0im, 0.0-1.0im]
roots(Poly([0, 0, 1]))     # [0.0, 0.0]
roots(poly([1,2,3,4]))     # [1.0,2.0,3.0,4.0]
```
"""
function roots{T}(p::Poly{T})
    R = promote_type(T, Float64)
    length(p) == 0 && return zeros(R, 0)

    num_leading_zeros = 0
    while abs(p[num_leading_zeros]) <= 2*eps(T)
        if num_leading_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_leading_zeros += 1
    end

    num_trailing_zeros = 0
    while abs(p[end - num_trailing_zeros]) <= 2*eps(T)
        num_trailing_zeros += 1
    end

    n = endof(p)-(num_leading_zeros + num_trailing_zeros)

    n < 1 && return zeros(R, length(p) - num_trailing_zeros - 1)

    companion = zeros(R, n, n)
    an = p[end-num_trailing_zeros]
    for i = 1:n-1
        companion[i,n] = -p[num_leading_zeros + i - 1] / an
        companion[i+1,i] = 1;
    end
    companion[end,end] = -p[end-num_trailing_zeros-1] / an
    D = eigvals(companion)
    r = zeros(eltype(D),length(p)-num_trailing_zeros-1)
    r[1:n] = D
    return r
end
roots{T}(p::Poly{Rational{T}}) = roots(convert(Poly{promote_type(T, Float64)}, p))

## compute gcd of two polynomials
"""

* `gcd(a::Poly, b::Poly)`: Finds the Greatest Common Denominator of
    two polynomials recursively using [Euclid's
    algorithm](http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm).

Example:
```
gcd(poly([1,1,2]), poly([1,2,3])) # returns (x-1)*(x-2)
```
"""
function gcd{T, S}(a::Poly{T}, b::Poly{S})
    if all(abs(b.a).<=2*eps(S))
        return a
    else
        s, r = divrem(a, b)
        return gcd(b, r)
    end
end


## Fit degree n polynomial to points
"""

`polyfit(x, y, n=length(x)-1, sym=:x )`: Fit a polynomial of degree
`n` through the points specified by `x` and `y` where `n <= length(x)
- 1` using least squares fit. When `n=length(x)-1` (the default), the
interpolating polynomial is returned. The optional fourth argument can
be used to pass the symbol for the returned polynomial.

Example:

```
xs = linspace(0, pi, 5)
ys = map(sin, xs)
polyfit(xs, ys, 2)
```

Original by [ggggggggg](https://github.com/Keno/Polynomials.jl/issues/19)
"""
function polyfit(x, y, n::Int=length(x)-1, sym::Symbol=:x)
    length(x) == length(y) || throw(DomainError)
    1 <= n <= length(x) - 1 || throw(DomainError)

    A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
    Poly(A \ y, sym)
end
polyfit(x,y,sym::Symbol) = polyfit(x,y,length(x)-1, sym)

### Pull in others
include("pade.jl")

end # module Poly
