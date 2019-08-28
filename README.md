# Polynomials

Basic arithmetic, integration, differentiation, evaluation, and root finding over dense univariate polynomials.

[![Polynomials](http://pkg.julialang.org/badges/Polynomials_0.6.svg)](http://pkg.julialang.org/?pkg=Polynomials)

Master branch:
[![Build Status](https://travis-ci.org/JuliaMath/Polynomials.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Polynomials.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaMath/Polynomials.jl/badge.svg)](https://coveralls.io/github/JuliaMath/Polynomials.jl)

Documentation:
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/Polynomials.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaMath.github.io/Polynomials.jl/latest)

#### Poly(a::Vector) where {T<:Number}

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

Construct a polynomial from its roots. This is in contrast to the
`Poly` constructor, which constructs a polynomial from its
coefficients.

```julia
// Represents (x-1)*(x-2)*(x-3)
julia> poly([1,2,3])
Poly(-6 + 11x - 6x^2 + x^3)
```

#### +, -, *, /, div, ==

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

julia> p * q
Poly(1 + 2x - x^2 - 2x^3)

julia> q / 2
Poly(0.5 - 0.5x^2)

julia> q ÷ p      # `div`, also `rem` and `divrem`
Poly(0.25 - 0.5x)
```

Note that operations involving polynomials with different variables will error.

```julia
julia> p = Poly([1, 2, 3], :x)
julia> q = Poly([1, 2, 3], :s)
julia> p + q
ERROR: Polynomials must have same variable.
```

To get the degree of the polynomial use the `degree` method

```
julia> degree(p)
2

julia> degree(p^2)
4

julia> degree(p-p)
-1
```

#### polyval(p::Poly, x::Number)

Evaluate the polynomial `p` at `x`.

```julia
julia> p = Poly([1, 0, -1])
julia> polyval(p, 0.1)
0.99
```

A call method is also available:

```julia
julia> p(0.1)
0.99
```


#### polyint(p::Poly, k::Number=0)

Integrate the polynomial `p` term by term, optionally adding constant
term `k`. The order of the resulting polynomial is one higher than the
order of `p`.

```julia
julia> polyint(Poly([1, 0, -1]))
Poly(x - 0.3333333333333333x^3)

julia> polyint(Poly([1, 0, -1]), 2)
Poly(2.0 + x - 0.3333333333333333x^3)
```

#### polyder(p::Poly)

Differentiate the polynomial `p` term by term. The order of the
resulting polynomial is one lower than the order of `p`.

```julia
julia> polyder(Poly([1, 3, -1]))
Poly(3 - 2x)
```

#### roots(p::Poly)

Return the roots (zeros) of `p`, with multiplicity. The number of
roots returned is equal to the order of `p`. By design, this is not type-stable,
the returned roots may be real or complex.

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

#### polyfit(x, y, n=length(x)-1)

* `polyfit`: fits a polynomial (of order `n`) to `x` and `y` using a least-squares approximation.

```julia
julia> xs = 1:4; ys = exp(xs); polyfit(xs, ys)
Poly(-7.717211620141281 + 17.9146616149694x - 9.77757245502143x^2 + 2.298404288652356x^3)
```

Visual example:

![newplot 42](https://user-images.githubusercontent.com/3156114/41799777-9ba00582-7627-11e8-94ef-15297ec8790e.png)

#### Other methods

Polynomial objects also have other methods:

* 0-based indexing is used to extract the coefficients of $a_0 + a_1
  x + a_2 x^2 + ...$, coefficients may be changed using indexing
  notation.

* `coeffs`: returns the entire coefficient vector

* `degree`: returns the polynomial degree, `length` is 1 plus the degree

* `variable`: returns the polynomial symbol as a degree 1 polynomial

* `norm`: find the `p`-norm of a polynomial

* `conj`: finds the conjugate of a polynomial over a complex fiel

* `truncate`: set to 0 all small terms in a polynomial; `chop` chops off
  any small leading values that may arise due to floating point
  operations.

* `gcd`: greatest common divisor of two polynomials.

* `Pade`: Return the
  [Pade approximant](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant)
  of order `m/n` for a polynomial as a `Pade` object.


## See also

* [MultiPoly.jl](https://github.com/daviddelaat/MultiPoly.jl) for sparse multivariate polynomials

* [MultivariatePolynomials.jl](https://github.com/blegat/MultivariatePolynomials.jl) for multivariate polynomials and moments of commutative or non-commutative variables

* [Nemo.jl](https://github.com/wbhart/Nemo.jl) for generic polynomial rings, matrix spaces, fraction fields, residue rings, power series

* [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl) for a fast complex polynomial root finder
