# Polynomial

Basic arithmetic, integration, differentiation, evaluation, and root finding over dense polynomials.

[![Build Status](https://travis-ci.org/vtjnash/Polynomial.jl.png?branch=master)](https://travis-ci.org/vtjnash/Polynomial.jl)

#### Poly{T<:Number}(a::Vector)
Construct a polynomial from its coefficients, highest order first.

```julia
julia> Poly([1,0,3,4])
Poly(1x^3 + 3x^1 + 4)
```

Leading zeros are stripped.

```julia
julia> Poly([0,1,2,3])
Poly(1x^2 + 2x^1 + 3)
```

#### poly(r::AbstractVector)
Construct a polynomial from its roots. This is in contrast to the `Poly` constructor, which constructs a polynomial from its coefficients.

```julia
// Represents (x - 1)*(x-2)*(x-3)
julia> poly([1,2,3])
Poly(1x^3 + -6x^2 + 11x^1 + -6)
```

#### +, -, *, /, ==

The usual arithmetic operators are overloaded to work on polynomials, and combinations of polynomials and scalars. Division by a polynomial is currently unimplemented. See [#1](https://github.com/vtjnash/Polynomial.jl/issues/1).

```julia
julia> a = Poly([1,2])
Poly(1x^1 + 2)

julia> b = Poly([1, 0, -1])
Poly(1x^2 + -1)

julia> 2*a
Poly(2x^1 + 4)

julia> 2 + a
Poly(1x^1 + 4)

julia> a - b
Poly(-1x^2 + 1x^1 + 3)

julia> a*b
Poly(1x^3 + 2x^2 + -1x^1 + -2)

julia> b/2
Poly(0.5x^2 + -0.5)

julia> a/b
ERROR: no method /(Poly{Int64}, Poly{Int64})
```

#### polyval(p::Poly, x::Number)
Evaluate the polynomial `p` at `x`.

```julia
julia> polyval(Poly([1, 0, -1]), 0.1)
-0.99
```

#### polyint(p::Poly, k::Number=0)
Integrate the polynomial `p` term by term, optionally adding constant term `k`. The order of the resulting polynomial is one higher than the order of `p`.

```julia
julia> polyint(Poly([1, 0, -1]))
Poly(0.3333333333333333x^3 + -1.0x^1)

julia> polyint(Poly([1, 0, -1]), 2)
Poly(0.3333333333333333x^3 + -1.0x^1 + 2.0)
```

#### polyder(p::Poly)
Differentiate the polynomial `p` term by term. The order of the resulting polynomial is one lower than the order of `p`.

```julia
julia> polyder(Poly([1, 3, -1]))
Poly(2x^1 + 3)
```

#### roots(p::Poly)
Return the roots (zeros) of `p`, with multiplicity. The number of roots returned is equal to the order of `p`. The returned roots may be real or complex.

```julia
julia> roots(Poly([1, 0, -1]))
2-element Array{Float64,1}:
 -1.0
  1.0

julia> roots(Poly([1, 0, 1]))
2-element Array{Complex{Float64},1}:
 0.0-1.0im
 0.0+1.0im

julia> roots(Poly([1, 0, 0]))
2-element Array{Float64,1}:
 0.0
 0.0
```

