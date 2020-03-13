# Extending Polynomials

The [`AbstractPolynomial`](@ref) type was made to be extended via a rich interface.

```@docs
AbstractPolynomial
```

A polynomial's  coefficients  are  relative to some *basis*. The `Polynomial` type relates coefficients  `[a0, a1,  ..., an]`, say,  to the  polynomial  `a0 +  a1*x + a2*x^  + ... +  an*x^n`,  through the standard  basis  `1,  x,  x^2, ..., x^n`.  New polynomial  types typically represent the polynomial through a different  basis. For example,  `CheyshevT` uses a basis  `T_0=1, T_1=x,  T_2=2x^2-1,  ...,  T_n  =  2xT_{n-1} - T_{n-2}`.  For this type  the  coefficients  `[a0,a1,...,an]` are associated with  the polynomial  `a0*T0  + a1*T_1 +  ...  +  an*T_n`.

To implement a new polynomial type, `P`, the following methods should
be implemented.

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
| `variable`| | Convenience to find monomial `x` in new  basis|

Check out both the [`Polynomial`](@ref) and [`ChebyshevT`](@ref) for examples of this interface being extended. [`Bernstein`](@ref) is an example where the basis depends on the  degree of  the polynomials being represented.
