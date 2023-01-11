# Chebyshev Polynomials

```@meta
DocTestSetup = quote
  using Polynomials
end
```


The [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials) are two sequences of polynomials, ``T_n`` and ``U_n``. The Chebyshev polynomials of the first kind, ``T_n``, can be defined by the recurrence relation:

```math
T_0(x)=1,\ T_1(x)=x
```

```math
T_{n+1}(x) = 2xT_n(x)-T_{n-1}(x)
```

The Chebyshev polynomioals of the second kind, ``U_n(x)``, can be defined by

```math
U_0(x)=1,\ U_1(x)=2x
```

```math
U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x)
```


Both ``T_n`` and ``U_n`` have degree ``n``, and any polynomial of degree ``n`` may be uniquely written as a linear combination of the polynomials ``T_0``, ``T_1``, ..., ``T_n`` (similarly with ``U_n``).


## First Kind

```@docs
ChebyshevT
```

The `ChebyshevT` type holds coefficients representing the polynomial ``a_0 T_0 + a_1 T_1 + ... + a_n T_n``.

For example, the basis polynomial ``T_4`` can be represented with `ChebyshevT([0,0,0,0,1])`.


### Conversion

[`ChebyshevT`](@ref) can be converted to [`Polynomial`](@ref) and vice-versa.

```jldoctest
julia> c = ChebyshevT([1, 0, 3, 4])
ChebyshevT(1⋅T_0(x) + 3⋅T_2(x) + 4⋅T_3(x))


julia> p = convert(Polynomial, c)
Polynomial(-2 - 12*x + 6*x^2 + 16*x^3)

julia> convert(ChebyshevT, p)
ChebyshevT(1.0⋅T_0(x) + 3.0⋅T_2(x) + 4.0⋅T_3(x))
```
