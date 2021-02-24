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

Check out both the [`Polynomial`](@ref) and [`ChebyshevT`](@ref) for examples of this interface being extended. 

The following shows a minimal example where the polynomial aliases the vector defining the coefficients. 
The constructor ensures that there are no trailing zeros. The method implemented below is the convenient  call syntax. This example uses the standard basis. For other bases, many more methods may be necessary to define  (again, refer to [`ChebyshevT`](@ref) for an example).

```jldoctest
julia> using Polynomials

julia> struct AliasPolynomial{T <: Number, X} <: Polynomials.StandardBasisPolynomial{T, X}
           coeffs::Vector{T}
           function AliasPolynomial{T, X}(coeffs::Vector{S}) where {T, X, S}
               N = findlast(!iszero, coeffs)
               N == nothing && return new{T,X}(zeros(T,1))
               resize!(coeffs, N)
               new{T,X}(coeffs)
           end
       end

julia> (p::CP{T})(x::S) where {T,S} = evalpoly(x, p.coeffs)

julia> Polynomials.@register AliasPolynomial
AliasPolynomial

julia> xs = [1,2,3,4];

julia> p = AliasPolynomial(xs)
AliasPolynomial(1 + 2*x + 3*x^2 + 4*x^3)

julia> q = AliasPolynomial(1.0, :y)
AliasPolynomial(1.0)

julia> p + q
AliasPolynomial(2.0 + 2.0*x + 3.0*x^2 + 4.0*x^3)

julia> p * p
AliasPolynomial(1 + 4*x + 10*x^2 + 20*x^3 + 25*x^4 + 24*x^5 + 16*x^6)

julia> (derivative ∘ integrate)(p) == p
true
```

For the `Polynomial` type, the default is to copy the array, so updates that don't change the size of the array. For this type, it might seem reasonable, to avoid allcoations, to update the coefficients in place for scalar addition and scalar multiplication. For that, the broadcast syntax might be useful. By default, the abstract polynomial type is not "broadcastable." To make this type allow scalare multiplation through the `.*=` syntax, we could do:

```jldoctest
julia> Base.broadcastable(p::AliasPolynomial) = p

julia> Base.ndims(::Type{<:AliasPolynomial}) = 1

julia> Base.copyto!(p::AliasPolynomial, x) = (copyto!(p.coeffs, x); chop!(p))
```

Then we might have:

```jldoctest
julia> p
AliasPolynomial(1 + 2*x + 3*x^2 + 4*x^3)

julia> p .*= 2
4-element Array{Int64,1}:
 2
 4
 6
 8

julia> p
AliasPolynomial(2 + 4*x + 6*x^2 + 8*x^3)
```

Though a bit cumbersome how the display of `p .+= 2` is the coefficients and not the polynomial, we can see that `p` is updated. 

However, scalar addition is defined differently than
