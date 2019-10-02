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

## Functions

```@docs
Poly
poly
degree
coeffs
variable
printpoly
polyval
polyint
polyder
polyfit
roots
```
