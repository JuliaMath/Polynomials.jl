# Polynomials.jl

Basic arithmetic, integration, differentiation, evaluation, and root finding over dense univariate polynomials.

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/Polynomials.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaMath.github.io/Polynomials.jl/dev)
[![Build Status](https://travis-ci.org/JuliaMath/Polynomials.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Polynomials.jl)
[![codecov](https://codecov.io/gh/JuliaMath/Polynomials.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/Polynomials.jl)


## Installation

```julia-repl
(v1.2) pkg> add Polynomials

julia> using Polynomials
```

## Usage

#### Available Polynomials

* `Polynomial` - Standard polynomials
* `ChebyshevT` - Chebyshev polynomials of the first kind

#### Construction and Evaluation

Construct a polynomial from its coefficients, lowest order first.

```julia-repl
julia> Polynomial([1,0,3,4])
Polynomial(1 + 3x^2 + 4x^3)
```

An optional variable parameter can be added.

```julia-repl
julia> Polynomial([1,2,3], :s)
Polynomial(1 + 2s + 3s^2)
```

Construct a polynomial from its roots.

```julia-repl
julia> fromroots([1,2,3]) # (x-1)*(x-2)*(x-3)
Polynomial(-6 + 11x - 6x^2 + x^3)
```

Evaluate the polynomial `p` at `x`.

```julia-repl
julia> p = Polynomial([1, 0, -1])
julia> p(0.1)
0.99
```

#### Arithmetic

The usual arithmetic operators are overloaded to work on polynomials, and combinations of polynomials and scalars.

```julia-repl
julia> p = Polynomial([1,2])
Polynomial(1 + 2x)

julia> q = Polynomial([1, 0, -1])
Polynomial(1 - x^2)

julia> 2p
Polynomial(2 + 4x)

julia> 2+p
Polynomial(3 + 2x)

julia> p - q
Poly(2x + x^2)

julia> p * q
Polynomial(1 + 2x - x^2 - 2x^3)

julia> q / 2
Polynomial(0.5 - 0.5x^2)

julia> q ÷ p  # `div`, also `rem` and `divrem`
Polynomial(0.25 - 0.5x)
```

Note that operations involving polynomials with different variables will error.

```julia-repl
julia> p = Polynomial([1, 2, 3], :x)
julia> q = Polynomial([1, 2, 3], :s)
julia> p + q
ERROR: Polynomials must have same variable.
```

#### Integrals and Derivatives

Integrate the polynomial `p` term by term, optionally adding constant
term `k`. The order of the resulting polynomial is one higher than the
order of `p`.

```julia-repl
julia> integral(Polynomial([1, 0, -1]))
Polynomial(x - 0.3333333333333333x^3)

julia> integral(Polynomial([1, 0, -1]), 2)
Polynomial(2.0 + x - 0.3333333333333333x^3)
```

Differentiate the polynomial `p` term by term. The order of the
resulting polynomial is one lower than the order of `p`.

```julia-repl
julia> derivative(Polynomial([1, 3, -1]))
Polynomial(3 - 2x)
```

#### Root-finding


Return the roots (zeros) of `p`, with multiplicity. The number of
roots returned is equal to the order of `p`. By design, this is not type-stable,
the returned roots may be real or complex.

```julia-repl
julia> roots(Polynomial([1, 0, -1]))
2-element Array{Float64,1}:
 -1.0
  1.0

julia> roots(Polynomial([1, 0, 1]))
2-element Array{Complex{Float64},1}:
 0.0+1.0im
 0.0-1.0im

julia> roots(Polynomial([0, 0, 1]))
2-element Array{Float64,1}:
 0.0
 0.0
```

#### Fitting arbitrary data

Fit a polynomial (of order `deg`) to `x` and `y` using a least-squares approximation.

```julia-repl
julia> xs = 1:4; ys = exp(xs);

julia> fit(xs, ys)
Polynomial(-7.717211620141281 + 17.9146616149694x - 9.77757245502143x^2 + 2.298404288652356x^3)
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

* `truncate`: set to 0 all small terms in a polynomial; 
* `chop` chops off any small leading values that may arise due to floating point operations.

* `gcd`: greatest common divisor of two polynomials.

* `Pade`: Return the
  [Pade approximant](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant) of order `m/n` for a polynomial as a `Pade` object.


## Related Packages

* [MultiPoly.jl](https://github.com/daviddelaat/MultiPoly.jl) for sparse multivariate polynomials

* [MultivariatePolynomials.jl](https://github.com/blegat/MultivariatePolynomials.jl) for multivariate polynomials and moments of commutative or non-commutative variables

* [Nemo.jl](https://github.com/wbhart/Nemo.jl) for generic polynomial rings, matrix spaces, fraction fields, residue rings, power series

* [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl) for a fast complex polynomial root finder


## Contributing

If you are interested in contributing, feel free to open an issue or pull request to get started.
