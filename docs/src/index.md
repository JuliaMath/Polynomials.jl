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
ERROR: ArgumentError: Polynomials have different indeterminates
[...]
```

Except for operations  involving constant polynomials.

```jldoctest
julia> p = Polynomial([1, 2, 3], :x)
Polynomial(1 + 2*x + 3*x^2)

julia> q = Polynomial(1, :y)
Polynomial(1)

julia> p+q
Polynomial(2 + 2*x + 3*x^2)
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
2-element Vector{Float64}:
 -1.0
  1.0

julia> roots(Polynomial([1, 0, 1]))
2-element Vector{ComplexF64}:
 0.0 - 1.0im
 0.0 + 1.0im

julia> roots(Polynomial([0, 0, 1]))
2-element Vector{Float64}:
 0.0
 0.0
```

For polynomials with multplicities, the non-exported `Polynomials.Multroot.multroot` function can avoid some numerical issues that `roots` will have. 

The `FactoredPolynomial` type stores the roots (with multiplicities) and the leading coefficent of a polynomial. In this example, the `multroot` function is used internally to identify the roots of `p` below, in the conversion from the `Polynomial` type to the `FactoredPolynomial` type:

```jldoctest
julia> p = Polynomial([24, -50, 35, -10, 1])
Polynomial(24 - 50*x + 35*x^2 - 10*x^3 + x^4)

julia> q = convert(FactoredPolynomial, p) # noisy form of `factor`:
FactoredPolynomial((x - 4.0000000000000036) * (x - 2.9999999999999942) * (x - 1.0000000000000002) * (x - 2.0000000000000018))
```


### Fitting arbitrary data

Fit a polynomial (of degree `deg`) to `x` and `y` using polynomial interpolation or a (weighted) least-squares approximation.

```@example
using Plots, Polynomials
xs = range(0, 10, length=10)
ys = @. exp(-xs)
f = fit(xs, ys) # degree = length(xs) - 1 
f2 = fit(xs, ys, 2) # degree = 2

scatter(xs, ys, markerstrokewidth=0, label="Data")
plot!(f, extrema(xs)..., label="Fit")
plot!(f2, extrema(xs)..., label="Quadratic Fit")
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
ChebyshevT(- 1.0⋅T_1(x) + 0.25⋅T_2(x) + 0.3333333333333333⋅T_3(x))

julia> convert(Polynomial, p1)
Polynomial(-2.0 + 2.0*x + 6.0*x^2)

julia> convert(ChebyshevT, Polynomial([1.0, 2,  3]))
ChebyshevT(2.5⋅T_0(x) + 2.0⋅T_1(x) + 1.5⋅T_2(x))
```

!!! warning
    The older  `Poly` type that this package used prior to `v0.7`  is implemented as an alternate basis  to provide support for older code bases. As of `v1.0`,  this type will be only available by executing `using Polynomials.PolyCompat`.

### Iteration

If its basis is implicit, then a polynomial may be  seen as just a vector of  coefficients. Vectors or 1-based, but, for convenience, polynomial types are 0-based, for purposes of indexing (e.g. `getindex`, `setindex!`, `eachindex`). Iteration over a polynomial steps through the underlying coefficients.

```jldoctest
julia> as = [1,2,3,4,5]; p = Polynomial(as);

julia> as[3], p[2], collect(p)[3]
(3, 3, 3)
```


The `pairs` iterator, iterates over the indices and coefficients, attempting to match how `pairs` applies to the underlying storage model:

```jldoctest
julia> v = [1,2,0,4]
4-element Vector{Int64}:
 1
 2
 0
 4

julia> p,ip,sp,lp = Polynomial(v), ImmutablePolynomial(v), SparsePolynomial(v), LaurentPolynomial(v, -1);

julia> collect(pairs(p))
4-element Vector{Pair{Int64, Int64}}:
 0 => 1
 1 => 2
 2 => 0
 3 => 4

julia> collect(pairs(ip)) == collect(pairs(p))
true

julia> collect(pairs(sp)) # unordered dictionary with only non-zero terms
3-element Vector{Pair{Int64, Int64}}:
 0 => 1
 3 => 4
 1 => 2

julia> collect(pairs(lp))
4-element Vector{Pair{Int64, Int64}}:
 -1 => 1
  0 => 2
  1 => 0
  2 => 4
```


The unexported `monomials` iterator iterates over the terms (`p[i]*Polynomials.basis(p,i)`) of the polynomial:

```jldoctest
julia> p = Polynomial([1,2,0,4], :u)
Polynomial(1 + 2*u + 4*u^3)

julia> collect(Polynomials.monomials(p))
4-element Vector{Any}:
 Polynomial(1)
 Polynomial(2*u)
 Polynomial(0)
 Polynomial(4*u^3)
```

The `map` function for polynomials is idiosyncratic, as iteration over
polynomials is over the vector of coefficients, but `map` will also
maintain the type of the polynomial. Here we use `map` to smooth out
the round-off error coming from the root-finding algorithm used
internally when converting to the `FactoredPolynomial` type:


```jldoctest
julia> p = Polynomial([24, -50, 35, -10, 1])
Polynomial(24 - 50*x + 35*x^2 - 10*x^3 + x^4)

julia> q = convert(FactoredPolynomial, p) # noisy form of `factor`:
FactoredPolynomial((x - 4.0000000000000036) * (x - 2.9999999999999942) * (x - 1.0000000000000002) * (x - 2.0000000000000018))

julia> map(round, q, digits=2)
FactoredPolynomial((x - 4.0) * (x - 2.0) * (x - 3.0) * (x - 1.0))
```


## Relationship between the `T` and `P{T,X}`

The addition of a polynomial and a scalar, such as

```jldoctest natural_inclusion
julia> using Polynomials

julia> p = Polynomial([1,2,3], :x)
Polynomial(1 + 2*x + 3*x^2)

julia> p + 3
Polynomial(4 + 2*x + 3*x^2)
```

seems natural, but in `Julia`, as `3` is of type `Int` and `p` of type `Polynomial{Int,:x}` some addition must be defined. The basic idea  is that `3` is promoted to the *constant* polynomial `3` with indeterminate `:x` as `3*one(p)` and then addition of `p + 3*one(p)` is performed.

This identification of a scalar with a constant polynomial can go both ways. If `q` is a *constant* polynomial of type `Polynomial{Int, :y}` then we should expect that `p+q` would be defined, as `p` plus the constant term of `q`. Indeed this is the case

```jldoctest natural_inclusion
julia> q = Polynomial(3, :y)
Polynomial(3)

julia> p + q
Polynomial(4 + 2*x + 3*x^2)
```

If `q` is non-constant, such as `variable(Polynomial, :y)`, then there would be an error due to the mismatched symbols. (The mathematical result would need a multivariate polynomial, not a univariate polynomial, as this package provides.)

The same conversion is done for polynomial multiplication: constant polynomials are treated as numbers; non-constant polynomials must have their symbols match. 

There is an oddity though the following two computations look the same, they are technically different:

```jldoctest natural_inclusion
julia> one(Polynomial, :x) + one(Polynomial, :y)
Polynomial(2.0)

julia> one(Polynomial, :y) + one(Polynomial, :x)
Polynomial(2.0)
```

Both are constant polynomials over `Int`, but the first has the indeterminate `:y`, the second `:x`. 

This technical difference causes no issues with polynomial addition or multiplication, as there constant polynomials are treated as numbers, but can be an issue when constant polynomials are used as array elements.

For arrays, the promotion of numbers to polynomials, allows natural constructions like:

```jldoctest natural_inclusion
julia> p = Polynomial([1,2],:x)
Polynomial(1 + 2*x)

julia> q = Polynomial([1,2],:y)  # non-constant polynomials with different indeterminates
Polynomial(1 + 2*y)

julia> [1 p]
1×2 Matrix{Polynomial{Int64, :x}}:
 Polynomial(1)  Polynomial(1 + 2*x)

julia> [1 one(q)]
1×2 Matrix{Polynomial{Int64, :y}}:
 Polynomial(1)  Polynomial(1)
```

However, as there would be an ambiguous outcome of the following

```jldoctest natural_inclusion
julia> [one(p) one(q)]
ERROR: ArgumentError: Polynomials have different indeterminates
[...]
```

an error is thrown.

In general, arrays with mixtures of non-constant polynomials with *different* indeterminates will error. By default, an error will occur when constant polynomials with different indeterminates are used as components. However, for *typed* arrays, conversion will allow such constructs to be used.

Using `one(q)` for a constant polynomial with indeterminate `:y` we have:

```jldoctest natural_inclusion
julia> P = typeof(p)
Polynomial{Int64, :x}

julia> P[one(p) one(q)]
1×2 Matrix{Polynomial{Int64, :x}}:
 Polynomial(1)  Polynomial(1)
```

Of course, by not being explicit, there are sill gotchas. For example, we can construct this matrix without a specific types:

```jldoctest natural_inclusion
julia> [one(p), one(q)+one(p)]
2-element Vector{Polynomial{Int64, :x}}:
 Polynomial(1)
 Polynomial(2)
```

but not this one:

```jldoctest natural_inclusion
julia> [one(p), one(p) + one(q)]
ERROR: ArgumentError: Polynomials have different indeterminates
[...]
```

Also, mixing types can result in unspecific symbols, as this example shows:

```jldoctest natural_inclusion
julia> [1 p; p 1] + [1 2one(q); 3 4] # array{P{T,:x}} + array{P{T,:y}}
2×2 Matrix{Polynomial{Int64, X} where X}:
 Polynomial(2)        Polynomial(3 + 2*x)
 Polynomial(4 + 2*x)  Polynomial(5)
```

Though were a non-constant polynomial with indeterminate `y` replacing
`2one(q)` above, that addition would throw an error.

## Rational functions

The package provides support for rational functions -- fractions of polynomials (for most types). The construction of the basic type mirrors the construction of rational numbers.

```jldoctest
julia> P = FactoredPolynomial
FactoredPolynomial

julia> p,q = fromroots(P, [1,2,3,4]), fromroots(P, [2,2,3,5])
(FactoredPolynomial((x - 4) * (x - 2) * (x - 3) * (x - 1)), FactoredPolynomial((x - 5) * (x - 2)² * (x - 3)))

julia> pq = p // q
((x - 4) * (x - 2) * (x - 3) * (x - 1)) // ((x - 5) * (x - 2)² * (x - 3))

julia> lowest_terms(pq)
((x - 4.0) * (x - 1.0)) // ((x - 5.0) * (x - 2.0))

julia> d,r = residues(pq); r
Dict{Float64, Vector{Float64}} with 2 entries:
  5.0 => [1.33333]
  2.0 => [0.666667]

julia> x = variable(p);

julia> for (λ, rs) ∈ r # reconstruct p/q from output of `residues`
           for (i,rᵢ) ∈ enumerate(rs)
               d += rᵢ//(x-λ)^i
           end
       end

julia> d
((x - 4.0) * (x - 1.0000000000000002)) // ((x - 5.0) * (x - 2.0))
```

A basic plot recipe is provided.

```@example
using Plots, Polynomials
P = FactoredPolynomial
p,q = fromroots(P, [1,2,3]), fromroots(P, [2,3,3,0])
plot(p//q)
savefig("rational_function.svg"); nothing # hide
```

![](rational_function.svg)


## Related Packages

* [StaticUnivariatePolynomials.jl](https://github.com/tkoolen/StaticUnivariatePolynomials.jl) Fixed-size univariate polynomials backed by a Tuple

* [MultiPoly.jl](https://github.com/daviddelaat/MultiPoly.jl) for sparse multivariate polynomials

* [DynamicPolynomals.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) Multivariate polynomials implementation of commutative and non-commutative variables

* [MultivariatePolynomials.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) for multivariate polynomials and moments of commutative or non-commutative variables

* [PolynomialRings.jl](https://github.com/tkluck/PolynomialRings.jl) A library for arithmetic and algebra with multi-variable polynomials.

* [AbstractAlgebra.jl](https://github.com/wbhart/AbstractAlgebra.jl), [Nemo.jl](https://github.com/wbhart/Nemo.jl) for generic polynomial rings, matrix spaces, fraction fields, residue rings, power series, [Hecke.jl](https://github.com/thofma/Hecke.jl) for algebraic number theory.

* [CommutativeAlgebra.jl](https://github.com/KlausC/CommutativeRings.jl) the start of a computer algebra system specialized to discrete calculations with support for polynomials.

* [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl) for a fast complex polynomial root finder. For larger degree problems, also [FastPolynomialRoots](https://github.com/andreasnoack/FastPolynomialRoots.jl) and [AMRVW](https://github.com/jverzani/AMRVW.jl).

* [SpecialPolynomials.jl](https://github.com/jverzani/SpecialPolynomials.jl) A package providing various polynomial types beyond the standard basis polynomials in `Polynomials.jl`. Includes interpolating polynomials, Bernstein polynomials, and classical orthogonal polynomials.

* [ClassicalOrthogonalPolynomials.jl](https://github.com/JuliaApproximation/ClassicalOrthogonalPolynomials.jl) A Julia package for classical orthogonal polynomials and expansions. Includes `chebyshevt`, `chebyshevu`, `legendrep`, `jacobip`, `ultrasphericalc`, `hermiteh`, and `laguerrel`. The same repository includes `FastGaussQuadrature.jl`, `FastTransforms.jl`, and the `ApproxFun` packages.



## Contributing

If you are interested in this project, feel free to open an issue or pull request! In general, any changes must be thoroughly tested, allow deprecation, and not deviate too far from the common interface. All PR's must have an updated project version, as well, to keep the continuous delivery cycle up-to-date.
