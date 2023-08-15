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

## Inspection

```@docs
degree
length
size
Polynomials.domain
mapdomain
chop
chop!
truncate
truncate!
iszero
Polynomials.isconstant
Polynomials.constantterm
isreal
real
isintegral
ismonic
Polynomials.hasnan
```

## Iteration

For the `Polynomial`  type, a natural mapping between the polynomial ``a_0 + a_1 x + a_2 x^2 + \cdots + a_n x^n`` with the coefficients ``(a_0, a_1, \dots, a_n))`` leads to the view point of a polynomial being a ``0``-based vector. Similarly, when the basis terms are not the standard basis. The `coeffs` method returns these coefficients in an iterable (a vector or tuple). For Laurent type polynomials, the coefficients between `firstindex(p)` and `lastindex(p)` are returned.

More generally, `pairs(p)` returns values `i => aᵢ` where the polynomial has terms ``a_i T_i`` for the basis ``T_i``. (For sparse polynomials these need not be in order and only terms where ``a_i \ne 0`` are given.) The `keys` and `values` methods iterate over `i` and `aᵢ`.

The `firstindex` method refers to the lowest stored basis index, which due to offsets need not be `0`. It will be no smaller than `Polynomials.minimumexponent`, which is the smalled allowed index for the polynomial type. The `lastindex` method refers to the last basis index. If the type allows trailing zeros (like `ImmutablePolynomial`) this will differ from the value returned by `degree`.

The `getindex(p,i)` method returns `p_i` or zero when out of bounds (if the element type of the polynomial has `zero(T)` defined).
For mutable polynomials, the `setindex!(p, val, i)` method sets `p[i]` to `val`. This may extend  the underlying storage container for some polynomial types. For `ImmutablePolynomial` the `@set!` macro from `Setfield` can be used with the typical `setindex!` notation.

The `map(fn, p)` method maps `fn` over the coefficients and returns a polynomial with the same polynomial type as `p`.

```@docs
coeffs
pairs
values
keys
firstindex
lastindex
eachindex
map
```

## Mathematical Functions

```@docs
zero
one
variable
Polynomials.basis
fromroots
gcd
derivative
integrate
roots
companion
fit
vander
```

## Plotting

Polynomials can be plotted directly using [Plots.jl](https://github.com/juliaplots/plots.jl) or [Makie.jl](https://github.com/MakieOrg/Makie.jl).

```julia
plot(::AbstractPolynomial; kwds...)
```

will automatically determine a range based on the critical points (roots, extrema and points of inflection).

```julia
plot(::AbstractPolynomial, a, b; kwds...)
```

will plot the polynomial within the range `[a, b]`.

### Example: The Polynomials.jl logo

```@example
using Plots, Polynomials
# T1, T2, T3, and T4:
chebs = [
  ChebyshevT([0, 1]),
  ChebyshevT([0, 0, 1]),
  ChebyshevT([0, 0, 0, 1]),
  ChebyshevT([0, 0, 0, 0, 1]),
]
colors = ["#4063D8", "#389826", "#CB3C33", "#9558B2"]

p = plot(legend=false, label="")
for (cheb, col) in zip(chebs, colors)
  plot!(cheb, c=col, lw=5)
end
savefig("chebs.svg"); nothing # hide
```

![](chebs.svg)
