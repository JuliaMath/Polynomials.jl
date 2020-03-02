# Bernstein Polynomials

```@meta
DocTestSetup = quote
  using Polynomials
end
```


The [Bernstein polynomials](https://en.wikipedia.org/wiki/Bernstein_polynomial) are a family of polynomials defined on the interval `[0,1]`. For each `n` there are `n+1` polynomials, given by: `b(n,nu) = choose(n,nu) x^nu * (1-x)^(n-nu)` for `nu` in `0:n`. Together, these form a basis that can represent any polynomial of degree `n` or less through a linear combination.


## Bernstein(N,T)

```@docs
Bernstein
Bernstein{N,T}()
```

The `Bernstein{N,T}` type holds coefficients representing the polynomial `a_0 b(n,0) + a_1 b(n,1) + ... + a_n b(n,n)`.

For example, using `n=3`, the monomial `x` is represented by `0 b(3,0) + 1/3 b(3,1) + 2/3 b(3,2) + 1  b(3,3)`, which can be constructed through `Bernstein{3, T}([0, 1/3, 2/3, 1])`


### Conversion

[`Bernsteing`](@ref) can be converted to [`Polynomial`](@ref) and vice-versa.

```jldoctest
julia> b = Bernstein([1, 0, 3, 4])
Bernstein(1⋅β(3, 0)(x) + 3⋅β(3, 2)(x) + 4⋅β(3, 3)(x))

julia> p = convert(Polynomial{Int}, b)
Polynomial(1 - 3*x + 12*x^2 - 6*x^3)

julia> convert(Bernstein{3, Float64},  p)
Bernstein(1.0⋅β(3, 0)(x) + 3.0⋅β(3, 2)(x) + 4.0⋅β(3, 3)(x))
```

Bernstein polynomials may be  converted  into a basis with a larger  degree:

```jldoctest
julia> convert(Bernstein{4, Float64}, b)
Bernstein(1.0⋅β(4, 0)(x) + 0.25⋅β(4, 1)(x) + 1.5⋅β(4, 2)(x) + 3.25⋅β(4, 3)(x) + 4.0⋅β(4, 4)(x))
```
