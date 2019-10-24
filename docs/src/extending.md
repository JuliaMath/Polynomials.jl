# Extending Polynomials

The [`AbstractPolynomial`](@ref) type was made to be extended via a rich interface. 

```@docs
AbstractPolynomial
```

To implement a new polynomial type, `P`, the following methods should be implemented. 

!!! note
    Promotion rules will always coerce towards the [`Polynomial`](@ref) type, so not all methods have to be implemented if you provide a conversion function.

As always, if the default implementation does not work or there are more efficient ways of implementing, feel free to overwrite functions from `common.jl` for your type.

| Function | Required | Notes |
|----------|:--------:|:------------|
| Constructor | x | |
| Type function (`(::P)(x)`) | x | |
| `convert(::Polynomial, ...)` | | Not required, but the library is built off the [`Polynomial`](@ref) type, so all operations are guaranteed to work with it. Also consider writing the inverse conversion method. |
| `domain` | x | Should return an  [`AbstractInterval`](https://invenia.github.io/Intervals.jl/stable/#Intervals-1) |
| `vander` | | Required for [`fit`](@ref) |
| `companion` | | Required for [`roots`](@ref) |
| `fromroots` | | By default, will form polynomials using `prod(variable(::P) - r)` for reach root `r`|
| `+(::P, ::P)` | | Addition of polynomials |
| `-(::P, ::P)` | | Subtraction of polynomials |
| `*(::P, ::P)` | | Multiplication of polynomials |
| `divrem` | | Required for [`gcd`](@ref)|

Check out both the [`Polynomial`](@ref) and [`ChebyshevT`](@ref) for examples of this interface being extended!
