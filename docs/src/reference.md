# Reference/API

All polynomials have the following functionality. In some cases, there is not a direct function call and therefore the polynomials have to be converted to the standard [`Polynomial`](@ref) type before continuing.

```@index
Pages = ["reference.md"]
```

```@meta
DocTestSetup = quote
  using Polynomials
end
```

## Inspection

```@docs
coeffs
order
degree
length
size
domain
chop
chop!
truncate
truncate!
```

## Arithmetic

All `AbstractPolynomials` have basic arithmetic operations defined on them (`+`, `-`, `*`, `/`, `÷`, `%`, `==`).

```jldoctest
julia> p = Polynomial([1, 2])
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
```

```@docs
gcd
```

## Mathematical Function

```@docs
fromroots
roots
derivative
integral
integrate
fit
companion
vander
```

## Plotting

Polynomials can be plotted directly using [Plots.jl](https://github.com/juliaplots/plots.jl).

```julia
plot(::AbstractPolynomial; kwds...)
```

will automatically determine a range based on the critical points (roots, extrema and points of inflection).

```julia
plot(::AbstractPolynomial, a, b; kwds...)
```

will plot the polynomial within the range `[a, b]`.
