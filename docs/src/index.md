# Polynomials.jl

Polynomials.jl is a Julia package that provides basic arithmetic, integration,
differentiation, evaluation, and root finding over dense univariate polynomials.

To install the package, run

```julia
Pkg.add("Polynomials")
```

The package can then be loaded into the current session using

```julia
using Polynomials
```

## Functions

```@meta
CurrentModule = Polynomials
```

```@docs
Poly
poly
degree
coeffs
variable
polyval
polyint
polyder
polyfit
roots
```
