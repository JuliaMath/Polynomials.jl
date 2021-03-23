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
| `Polynomials.evalpoly(x, p::P)` |  to evaluate the polynomial at `x` (`Base.evalpoly` okay post `v"1.4.0"`) |
| `domain` | x | Should return an  [`AbstractInterval`](https://invenia.github.io/Intervals.jl/stable/#Intervals-1) |
| `vander` | | Required for [`fit`](@ref) |
| `companion` | | Required for [`roots`](@ref) |
| `*(::P, ::P)` | | Multiplication of polynomials |
| `divrem` | | Required for [`gcd`](@ref)|
| `one`| | Convenience to find constant in new basis |
| `variable`| | Convenience to find monomial `x` in new  basis|

Check out both the [`Polynomial`](@ref) and [`ChebyshevT`](@ref) for examples of this interface being extended. 

## Example

The following shows a minimal example where the polynomial aliases the vector defining the coefficients. 
The constructor ensures that there are no trailing zeros. The `@register` call ensures a common interface. This example subtypes `StandardBasisPolynomial`, not `AbstractPolynomial`, and consequently inherits the methods above that otherwise would have been required. For other bases, more methods may be necessary to define (again, refer to [`ChebyshevT`](@ref) for an example).

```jldoctest AliasPolynomial
julia> using Polynomials

julia> struct AliasPolynomial{T <: Number, X} <: Polynomials.StandardBasisPolynomial{T, X}
                  coeffs::Vector{T}
                  function AliasPolynomial{T, X}(coeffs::Vector{S}) where {T, X, S}
                      p = new{T,X}(coeffs)
                      chop!(p)
                  end
              end

julia> Polynomials.@register AliasPolynomial
```

To see this new polynomial type in action, we have:

```jldoctest AliasPolynomial
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

julia> p(3)
142
```

For the `Polynomial` type, the default on operations is to copy the array. For this type, it might seem reasonable -- to avoid allocations -- to update the coefficients in place for scalar addition and scalar multiplication. 

Scalar addition, `p+c`, defaults to `p + c*one(p)`, or polynomial addition, which is not inplace without addition work. As such, we create a new method and an infix operator

```jldoctest AliasPolynomial
julia> function scalar_add!(p::AliasPolynomial{T}, c::T) where {T}
           p.coeffs[1] += c
           p
       end;

julia> p::AliasPolynomial ⊕ c::Number = scalar_add!(p,c);

```

Then we have:

```jldoctest AliasPolynomial
julia> p
AliasPolynomial(1 + 2*x + 3*x^2 + 4*x^3)

julia> p ⊕ 2
AliasPolynomial(3 + 2*x + 3*x^2 + 4*x^3)

julia> p
AliasPolynomial(3 + 2*x + 3*x^2 + 4*x^3)
```

The viewpoint that a polynomial represents a vector of coefficients  leads to an expectation that vector operations should match when possible. Scalar multiplication is a vector operation, so it seems reasonable to override the broadcast machinery to implement an in place operation (e.g. `p .*= 2`). By default, the polynomial types are not broadcastable over their coefficients. We would need to make a change there and modify the `copyto!` function:


```jldoctest AliasPolynomial
julia> Base.broadcastable(p::AliasPolynomial) = p.coeffs;


julia> Base.ndims(::Type{<:AliasPolynomial}) = 1


julia> Base.copyto!(p::AliasPolynomial, x) = (copyto!(p.coeffs, x); chop!(p));

```

The last `chop!` call would ensure that there are no trailing zeros in the coefficient vector after multiplication, as multiplication by `0` is possible.

Then we might have:

```jldoctest AliasPolynomial
julia> p
AliasPolynomial(3 + 2*x + 3*x^2 + 4*x^3)

julia> p .*= 2
AliasPolynomial(6 + 4*x + 6*x^2 + 8*x^3)

julia> p
AliasPolynomial(6 + 4*x + 6*x^2 + 8*x^3)

julia> p ./= 2
AliasPolynomial(3 + 2*x + 3*x^2 + 4*x^3)
```

Trying to divide again would throw an error, as the result would not fit with the integer type of `p`. 

Now `p` is treated as the vector `p.coeffs`, as regards broadcasting, so some things may be surprising, for example this expression returns a vector, not a polynomial:

```jldoctest AliasPolynomial
julia> p .+ 2
4-element Array{Int64,1}:
 5
 4
 5
 6
```

