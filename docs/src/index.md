# Polynomials.jl

Polynomials.jl is a Julia package that provides basic arithmetic, integration,
differentiation, evaluation, and root finding over dense univariate polynomials.

To install the package, run

```julia
(v1.2) pkg> add Polynomials
```

The package can then be loaded into the current session using

```julia
using Polynomials
```

```@meta
DocTestSetup = quote
  using Polynomials
end
```

## Quick Start

### Construction and Evaluation

Construct a polynomial from its coefficients, lowest order first.

```jldoctest
julia> Polynomial([1,0,3,4])
Polynomial(1 + 3*x^2 + 4*x^3)
```

An optional variable parameter can be added.

```jldoctest
julia> Polynomial([1,2,3], :s)
Polynomial(1 + 2*s + 3*s^2)
```

Construct a polynomial from its roots.

```jldoctest
julia> fromroots([1,2,3]) # (x-1)*(x-2)*(x-3)
Polynomial(-6 + 11*x - 6*x^2 + x^3)
```

Evaluate the polynomial `p` at `x`.

```jldoctest
julia> p = Polynomial([1, 0, -1])
Polynomial(1 - x^2)

julia> p(1)
0

```

### Arithmetic

The usual arithmetic operators are overloaded to work on polynomials, and combinations of polynomials and scalars.

```jldoctest
julia> p = Polynomial([1,2])
Polynomial(1 + 2*x)

julia> q = Polynomial([1, 0, -1])
Polynomial(1 - x^2)

julia> 2p
Polynomial(2 + 4*x)

julia> 2 + p
Polynomial(3 + 2*x)

julia> p - q
Polynomial(2*x + x^2)

julia> p * q
Polynomial(1 + 2*x - x^2 - 2*x^3)

julia> q / 2
Polynomial(0.5 - 0.5*x^2)

julia> q ÷ p  # `div`, also `rem` and `divrem`
Polynomial(0.25 - 0.5*x)
```

Note that operations involving polynomials with different variables will error.

```jldoctest
julia> p = Polynomial([1, 2, 3], :x)
Polynomial(1 + 2*x + 3*x^2)

julia> q = Polynomial([1, 2, 3], :s)
Polynomial(1 + 2*s + 3*s^2)

julia> p + q
ERROR: Polynomials must have same variable
[...]
```

### Integrals and Derivatives

Integrate the polynomial `p` term by term, optionally adding constant
term `C`. The degree of the resulting polynomial is one higher than the
degree of `p`.

```jldoctest
julia> integrate(Polynomial([1, 0, -1]))
Polynomial(1.0*x - 0.3333333333333333*x^3)

julia> integrate(Polynomial([1, 0, -1]), 2)
Polynomial(2.0 + 1.0*x - 0.3333333333333333*x^3)
```

Differentiate the polynomial `p` term by term. The degree of the
resulting polynomial is one lower than the degree of `p`.

```jldoctest
julia> derivative(Polynomial([1, 3, -1]))
Polynomial(3 - 2*x)
```

### Root-finding

Return the roots (zeros) of `p`, with multiplicity. The number of
roots returned is equal to the order of `p`. By design, this is not type-stable,
the returned roots may be real or complex.

```jldoctest
julia> roots(Polynomial([1, 0, -1]))
2-element Array{Float64,1}:
  1.0
 -1.0

julia> roots(Polynomial([1, 0, 1]))
2-element Array{Complex{Float64},1}:
  0.0 - 1.0im
  0.0 + 1.0im

julia> roots(Polynomial([0, 0, 1]))
2-element Array{Float64,1}:
  0.0
 -0.0
```

### Fitting arbitrary data

Fit a polynomial (of degree `deg`) to `x` and `y` using polynomial interpolation or a (weighted) least-squares approximation.

```@example
using Plots, Polynomials
xs = range(0, 10, length=10)
ys = exp.(-xs)
f = fit(xs, ys)  # fit(xs, ys, k)  for fitting a kth degreee polynomial

scatter(xs, ys, label="Data");
plot!(f, extrema(xs)..., label="Fit");
savefig("polyfit.svg"); nothing # hide
```

![](polyfit.svg)

### Other bases

A polynomial, e.g. `a_0 + a_1 x + a_2 x^2 + ... + a_n x^n`, can be seen as a collection of coefficients, `[a_0, a_1, ..., a_n]`, relative to some polynomial basis. The most  familiar basis being  the standard one: `1`, `x`, `x^2`, ...  Alternative bases are possible.  The `ChebyshevT` polynomials are  implemented, as an example. Instead of `Polynomial`  or  `Polynomial{T}`, `ChebyshevT` or  `ChebyshevT{T}` constructors are used:

```jldoctest
julia> p1 = ChebyshevT([1.0, 2.0, 3.0])
ChebyshevT(1.0⋅T_0(x) + 2.0⋅T_1(x) + 3.0⋅T_2(x))

julia> p2 = ChebyshevT{Float64}([0, 1, 2])
ChebyshevT(1.0⋅T_1(x) + 2.0⋅T_2(x))

julia> p1 + p2
ChebyshevT(1.0⋅T_0(x) + 3.0⋅T_1(x) + 5.0⋅T_2(x))

julia> p1 * p2
ChebyshevT(4.0⋅T_0(x) + 4.5⋅T_1(x) + 3.0⋅T_2(x) + 3.5⋅T_3(x) + 3.0⋅T_4(x))

julia> derivative(p1)
ChebyshevT(2.0⋅T_0(x) + 12.0⋅T_1(x))

julia> integrate(p2)
ChebyshevT(0.25⋅T_0(x) - 1.0⋅T_1(x) + 0.25⋅T_2(x) + 0.3333333333333333⋅T_3(x))

julia> convert(Polynomial, p1)
Polynomial(-2.0 + 2.0*x + 6.0*x^2)

julia> convert(ChebyshevT, Polynomial([1.0, 2,  3]))
ChebyshevT(2.5⋅T_0(x) + 2.0⋅T_1(x) + 1.5⋅T_2(x))
```


### Iteration

If its basis is implicit, then a polynomial may be  seen as just a vector of  coefficients. Vectors or 1-based, but, for convenience, polynomial types are 0-based, for purposes of indexing (e.g. `getindex`, `setindex!`, `eachindex`). Iteration over a polynomial steps through the basis vectors, e.g. `a_0`, `a_1*x`, ...

```jldoctest
julia> as = [1,2,3,4,5]; p = Polynomial(as);

julia> as[3], p[2], collect(p)[3]
(3, 3, Polynomial(3*x^2))
```

## Related Packages

* [MultiPoly.jl](https://github.com/daviddelaat/MultiPoly.jl) for sparse multivariate polynomials

* [MultivariatePolynomials.jl](https://github.com/blegat/MultivariatePolynomials.jl) for multivariate polynomials and moments of commutative or non-commutative variables

* [AbstractAlgeebra.jl](https://github.com/wbhart/AbstractAlgebra.jl) for generic polynomial rings, matrix spaces, fraction fields, residue rings, power series.

* [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl) for a fast complex polynomial root finder

## Contributing

If you are interested in this project, feel free to open an issue or pull request! In general, any changes must be thoroughly tested, allow deprecation, and not deviate too far from the common interface. All PR's must have an updated project version, as well, to keep the continuous delivery cycle up-to-date.
