# Extending Polynomials

The [`AbstractPolynomial`](@ref) type was made to be extended via a rich interface. To implement a new polynomial type, the following methods should be implemented. 

```@docs
AbstractPolynomial
```

| Function | Required | Description |
|----------|:--------:|:------------|
| Constructor | x | |
| Type function (`(::MyType)()`) | x | | 
| `convert(::Polynomial, ...)` | | Not required, but the library is built off the [`Polynomial`](@ref) type, so all operations are guaranteed to work with it. Also consider writing the inverse conversion method. |
| `domain` | x | |
| `vander` | | Required for [`fit`](@ref) |
| `companion` | | Required for [`roots`](@ref) |
| `fromroots` | | |
| `+` | | | 
| `-` | | |
| `*` | | |
| `divrem` | | |


