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

## Quick Start

```jldoctest
julia> p = fromroots([-1, 3])
Polynomial(-3 - 2*x + x^2)

julia> p.(-1:1)
3-element Array{Int64,1}:
  0
 -3
 -4

julia> roots(p)
2-element Array{Float64,1}:
  3.0000000000000004
 -0.9999999999999998

julia> pint = integral(p)
Polynomial(-3.0*x - 1.0*x^2 + 0.3333333333333333*x^3)

julia> pder = derivative(pint)
Polynomial(-3.0 - 2.0*x + 1.0*x^2)

julia> integrate(p, 0, 1) == pint(1) - pint(0)
true

julia> x = variable()
Polynomial(1.0*x)

julia> p2 = 100 + 2x - 3x^2 + x^3
Polynomial(100.0 + 2.0*x - 3.0*x^2 + 1.0*x^3)

julia> p + p2
Polynomial(97.0 - 2.0*x^2 + 1.0*x^3)

julia> p * x
Polynomial(-300.0 - 206.0*x + 105.0*x^2 + 5.0*x^3 - 5.0*x^4 + 1.0*x^5)

julia> p2 ÷ p
Polynomial(-1.0 + 1.0*x)

julia> p2 % p
Polynomial(97.0 + 3.0*x)
```

See [Usage](@ref) for more details.

## Extending

If you want to implement your own polynomial, check out the guide in [Extending Polynomials](@ref) to see how to get started.

## Contributing

If you are interested in this project, feel free to open an issue or pull request! In general, any changes must be thoroughly tested, allow deprecation, and not deviate too far from the common interface. All PR's must have an updated project version, as well, to keep the continuous delivery cycle up-to-date.
