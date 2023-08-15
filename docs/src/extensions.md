# Extensions

As of `v1.9` of `Julia`, packages can provide extension code which is loaded when external packages are loaded.

## Makie

When `Makie` is loaded, a plot recipe is provided.

## ChainRulesCore

When `ChainRulesCore` is loaded, a `frule` and `rrule` is defined for to integrate with different autodifferentiation packages.

## MutableArithmetics

When the `MutableArithmetics` package is loaded, an extension provides its functionality for a few polynomial types, described in the following. Prior to `v1.9` the external package `PolynomialsMutableArithmetics` provided the same functionality.


While polynomials of type `Polynomial` are mutable objects, operations such as
`+`, `-`, `*`, always create new polynomials without modifying its arguments.
The time needed for these allocations and copies of the polynomial coefficients
may be noticeable in some use cases. This is amplified when the coefficients
are for instance `BigInt` or `BigFloat` which are mutable themselves.
This can be avoided by modifying existing polynomials to contain the result
of the operation using the [MutableArithmetics (MA) API](https://github.com/jump-dev/MutableArithmetics.jl).

Consider for instance the following arrays of polynomials

```julia
using Polynomials
d, m, n = 30, 20, 20
p(d) = Polynomial(big.(1:d))
A = [p(d) for i in 1:m, j in 1:n]
b = [p(d) for i in 1:n]
```

In this case, the arrays are mutable objects for which the elements are mutable
polynomials which have mutable coefficients (`BigInt`s).
These three nested levels of mutable objects communicate with the MA
API in order to reduce allocation.
Calling `A * b` requires approximately 40 MiB due to 2 M allocations
as it does not exploit any mutability.

```julia
using MutableArithmetics  # or `using PolynomialsMutableArithmetics` to register `Polynomials` with `MutableArithmetics`

const MA = MutableArithmetics
MA.operate(*, A, b)
```

exploits the mutability and hence only allocates approximately 70 KiB due to 4 k
allocations.

If the resulting vector is already allocated, e.g.,

```julia
z(d) = Polynomial([zero(BigInt) for i in 1:d])
c = [z(2d - 1) for i in 1:m]
```

then we can exploit its mutability with

```julia
MA.operate!(MA.add_mul, c, A, b)
```

to reduce the allocation down to 48 bytes due to 3 allocations.

These remaining allocations are due to the `BigInt` buffer used to
store the result of intermediate multiplications. This buffer can be
preallocated with:

```julia
buffer = MA.buffer_for(MA.add_mul, typeof(c), typeof(A), typeof(b))
MA.buffered_operate!(buffer, MA.add_mul, c, A, b)
```

then the second line is allocation-free.

The `MA.@rewrite` macro rewrite an expression into an equivalent code that
exploit the mutability of the intermediate results.
For instance
```julia
MA.@rewrite(A1 * b1 + A2 * b2)
```
is rewritten into
```julia
c = MA.operate!(MA.add_mul, MA.Zero(), A1, b1)
MA.operate!(MA.add_mul, c, A2, b2)
```
which is equivalent to
```julia
c = MA.operate(*, A1, b1)
MA.mutable_operate!(MA.add_mul, c, A2, b2)
```

!!! note
    Note that currently, only the `Polynomial` and `Polynomials.PnPolynomial` types implement the API and  only
part of it is implemented

## PolyCompat

While not an extension, the older  `Poly` type that this package used prior to `v0.7`  is implemented as an alternate basis
and provided on an opt-in bases by executing `using Polynomials.PolyCompat`. This is to provide support for older code bases.
