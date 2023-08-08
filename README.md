# Polynomials.jl

Basic arithmetic, integration, differentiation, evaluation, root finding, and fitting for
univariate [polynomials](https://en.wikipedia.org/wiki/Polynomial).

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/Polynomials.jl/stable)
[![CI](https://github.com/JuliaMath/Polynomials.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaMath/Polynomials.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaMath/Polynomials.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/Polynomials.jl)



## Available Types of Polynomials

There are a number of different standard basis polynomials primarily distinguished by their storage type:

* `Polynomial` –⁠ standard basis polynomials, $a(x) = a_0 + a_1 x + a_2 x^2 + … + a_n  x^n$ for $n ≥ 0$.
* `ImmutablePolynomial` –⁠ standard basis polynomials backed by a [Tuple type](https://docs.julialang.org/en/v1/manual/functions/#Tuples-1) for faster evaluation of values
* `SparsePolynomial` –⁠ standard basis polynomial backed by a [dictionary](https://docs.julialang.org/en/v1/base/collections/#Dictionaries-1) to hold  sparse high-degree  polynomials
* `LaurentPolynomial` –⁠ [Laurent polynomials](https://docs.julialang.org/en/v1/base/collections/#Dictionaries-1), $a(x) = a_m x^m + … + a_n x^n$ for $m ≤ n$ and $m,n ∈ ℤ$. This is backed by an [offset array](https://github.com/JuliaArrays/OffsetArrays.jl); for example, if $m<0$ and $n>0$, we obtain $a(x) = a_m x^m + … + a_{-1} x^{-1} + a_0 + a_1 x + … +  a_n x^n$
* `FactoredPolynomial` –⁠ standard basis polynomials, storing the roots, with multiplicity, and leading coefficient of a polynomial

The  `ChebyshevT` –⁠ [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials) of the first kind illustrate how other bases may be employed with the package.

In addition `RationalFunction` is a type for ratios of polynomials.


### Fitting arbitrary data

Fit a polynomial (of degree `deg` or less) to `x` and `y` using a least-squares approximation.

```julia
julia> xs = 0:4; ys = @. exp(-xs) + sin(xs);

julia> fit(xs, ys) |> p -> round.(coeffs(p), digits=4) |> Polynomial
Polynomial(1.0 + 0.0593*x + 0.3959*x^2 - 0.2846*x^3 + 0.0387*x^4)

julia> fit(ChebyshevT, xs, ys, 2) |> p -> round.(coeffs(p), digits=4) |> ChebyshevT
ChebyshevT(0.5413⋅T_0(x) - 0.8991⋅T_1(x) - 0.4238⋅T_2(x))
```

Visual example:

![fit example](https://user-images.githubusercontent.com/14099459/70382587-9e055500-1902-11ea-8952-3f03ae08b7dc.png)

### Other methods

Polynomial objects also have other methods:

* For standard basis polynomials, 0-based indexing is used to extract
  the coefficients of `[a0, a1, a2, ...]`; for mutable polynomials,
  coefficients may be changed using indexing notation.

* `coeffs`: returns the coefficients

* `degree`: returns the polynomial degree, `length` is number of stored coefficients

* `variable`: returns the polynomial symbol as a polynomial in the underlying type

* `LinearAlgebra.norm`: find the `p`-norm of a polynomial

* `conj`: finds the conjugate of a polynomial over a complex field

* `truncate`: set to 0 all small terms in a polynomial;

* `chop` chops off any small leading values that may arise due to floating point operations.

* `gcd`: greatest common divisor of two polynomials.

* `Pade`: Return the
  [Padé approximant](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant) of order `m/n` for a polynomial as a `Pade` object.


## Related Packages

* [StaticUnivariatePolynomials.jl](https://github.com/tkoolen/StaticUnivariatePolynomials.jl) Fixed-size univariate polynomials backed by a Tuple

* [MultiPoly.jl](https://github.com/daviddelaat/MultiPoly.jl) for sparse multivariate polynomials

* [DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) Multivariate polynomials implementation of commutative and non-commutative variables

* [MultivariatePolynomials.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) for multivariate polynomials and moments of commutative or non-commutative variables

* [PolynomialRings.jl](https://github.com/tkluck/PolynomialRings.jl) A library for arithmetic and algebra with multi-variable polynomials.

* [AbstractAlgebra.jl](https://github.com/wbhart/AbstractAlgebra.jl), [Nemo.jl](https://github.com/wbhart/Nemo.jl) for generic polynomial rings, matrix spaces, fraction fields, residue rings, power series, [Hecke.jl](https://github.com/thofma/Hecke.jl) for algebraic number theory.

* [LaurentPolynomials.jl](https://github.com/jmichel7/LaurentPolynomials.jl) A package for Laurent polynomials.

* [CommutativeAlgebra.jl](https://github.com/KlausC/CommutativeRings.jl) the start of a computer algebra system specialized to discrete calculations with support for polynomials.

* [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl) for a fast complex polynomial root finder. For larger degree problems, also [FastPolynomialRoots](https://github.com/andreasnoack/FastPolynomialRoots.jl) and [AMRVW](https://github.com/jverzani/AMRVW.jl). For real roots only [RealPolynomialRoots](https://github.com/jverzani/RealPolynomialRoots.jl).


* [SpecialPolynomials.jl](https://github.com/jverzani/SpecialPolynomials.jl) A package providing various polynomial types beyond the standard basis polynomials in `Polynomials.jl`. Includes interpolating polynomials, Bernstein polynomials, and classical orthogonal polynomials.

* [ClassicalOrthogonalPolynomials.jl](https://github.com/JuliaApproximation/ClassicalOrthogonalPolynomials.jl) A Julia package for classical orthogonal polynomials and expansions. Includes `chebyshevt`, `chebyshevu`, `legendrep`, `jacobip`, `ultrasphericalc`, `hermiteh`, and `laguerrel`. The same repository includes `FastGaussQuadrature.jl`, `FastTransforms.jl`, and the `ApproxFun` packages.


## Legacy code

As of v0.7, the internals of this package were greatly generalized and new types and method names were introduced. For compatibility purposes, legacy code can be run after issuing `using Polynomials.PolyCompat`.

## Contributing

If you are interested in contributing, feel free to open an issue or pull request to get started.
