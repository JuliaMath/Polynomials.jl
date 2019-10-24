# Chebyshev Polynomials

```@meta
DocTestSetup = quote
  using Polynomials
end
```

## First Kind

```@docs
ChebyshevT
ChebyshevT()
```

### Conversion

[`ChebyshevT`](@ref) can be converted to [`Polynomial`](@ref) and vice-versa.

```jldoctest
julia> c = ChebyshevT([1, 0, 3, 4])
ChebyshevT([1, 0, 3, 4])

julia> p = convert(Polynomial, c)
Polynomial(-2 - 12*x + 6*x^2 + 16*x^3)

julia> convert(ChebyshevT, p)
ChebyshevT([1.0, 0.0, 3.0, 4.0])
```
