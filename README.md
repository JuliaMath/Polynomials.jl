# Polynomials

Basic arithmetic, integration, differentiation, evaluation, and root finding over dense univariate polynomials.

[![Build Status](https://travis-ci.org/Keno/Polynomials.jl.png?branch=master)](https://travis-ci.org/Keno/Polynomials.jl)

#### Poly{T<:Number}(a::Vector)
Construct a polynomial from its coefficients, lowest order first.
```julia
julia> Poly([1,0,3,4])
Poly(1 + 3x^2 + 4x^3)
```

An optional variable parameter can be added.
```julia
julia> Poly([1,2,3], :s)
Poly(1 + 2s + 3s^2)
```

#### poly(r::AbstractVector)
Construct a polynomial from its roots. This is in contrast to the `Poly` constructor, which constructs a polynomial from its coefficients.
```julia
// Represents (x-1)*(x-2)*(x-3)
julia> poly([1,2,3])
Poly(-6 + 11x - 6x^2 + x^3)
```

#### +, -, *, /, ==

The usual arithmetic operators are overloaded to work on polynomials, and combinations of polynomials and scalars.
```julia
julia> p = Poly([1,2])
Poly(1 + 2x)

julia> q = Poly([1, 0, -1])
Poly(1 - x^2)

julia> 2p
Poly(2 + 4x)

julia> 2+p
Poly(3 + 2x)

julia> p - q
Poly(2x + x^2)

julia> p*q
Poly(1 + 2x - x^2 - 2x^3)

julia> q/2
Poly(0.5 - 0.5x^2)
```

Note that operations involving polynomials with different variables will error.

To get the degree of the polynomial use `degree` method
```
julia> degree(p)
1

julia> degree(p^2)
2

julia> degree(p-p)
0
```

```julia
julia> p = Poly([1, 2, 3], :x)
julia> q = Poly([1, 2, 3], :s)
julia> p + q
ERROR: Polynomials must have same variable.
```

#### polyval(p::Poly, x::Number)
Evaluate the polynomial `p` at `x`.

```julia
julia> polyval(Poly([1, 0, -1]), 0.1)
0.99
```

#### polyint(p::Poly, k::Number=0)
Integrate the polynomial `p` term by term, optionally adding constant term `k`. The order of the resulting polynomial is one higher than the order of `p`.
```julia
julia> polyint(Poly([1, 0, -1]))
Poly(x - 0.3333333333333333x^3)

julia> polyint(Poly([1, 0, -1]), 2)
Poly(2.0 + x - 0.3333333333333333x^3)
```

#### polyder(p::Poly)
Differentiate the polynomial `p` term by term. The order of the resulting polynomial is one lower than the order of `p`.
```julia
julia> polyder(Poly([1, 3, -1]))
Poly(3 - 2x)
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
 0.0+1.0im
 0.0-1.0im

julia> roots(Poly([0, 0, 1]))
2-element Array{Float64,1}:
 0.0
 0.0
```
