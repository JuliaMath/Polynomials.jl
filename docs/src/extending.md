# Extending Polynomials

The [`AbstractUnivariatePolynomial`](@ref) type was made to be extended.

A polynomial's  coefficients  are  relative to some *basis*. The `Polynomial` type relates coefficients  `[a0, a1,  ..., an]`, say,  to the  polynomial  ``a_0 +  a_1\cdot x + a_2\cdot x^2  + \cdots +  a_n\cdot x^n``,  through the standard  basis  ``1,  x,  x^2, ..., x^n``.  New polynomial  types typically represent the polynomial through a different  basis. For example,  `CheyshevT` uses a basis  ``T_0=1, T_1=x,  T_2=2x^2-1,  \cdots,  T_n  =  2xT_{n-1} - T_{n-2}``.  For this type  the  coefficients  `[a0,a1,...,an]` are associated with  the polynomial  ``a0\cdot T_0  + a_1 \cdot T_1 +  \cdots  +  a_n\cdot T_n`.

A polynomial type consists of a container type (with parent type `AbstractUnivariatePolynomial`) and a basis type (with parent type `AbstractBasis`). There a several different storage types implemented.

To implement a new polynomial type, `P`, the following methods should be implemented:


| Function | Required | Notes |
|----------|:--------:|:------------|
| A container type           | x | Usually selected from an available one. |
| A basis type               | x | |
| `variable`                 |   | Convenience to find the monomial `x` in the new basis.|
| `Base.evalpoly(x, p::P)`   | x | To evaluate the polynomial at `x` |
| `*(::P, ::P)`              |   | Multiplication of polynomials |
| `convert(::P, p::Polynomial)`  | Defaults to polynomial evaluation. Can be used to define `*` by round trip through `Polynomial` type|
| `convert(::Polynomial, p)` |   | Defaults to polynomial evaluation, which uses `evalpoly`, `variable`, `*`|
| `scalar_add(c::S, ::P)`    |   | Scalar addition. Default requires `one` to be defined. |
| `one`                      | x | Convenience to find the constant $1$ in the new basis.  |
| `map(f, p)`                | x | Used to define scalar multiplication |
| `divrem`                   |   | Required for [`gcd`](@ref)|
| `vander`                   |   | Required for [`fit`](@ref) |
| `companion`                |   | Required for [`roots`](@ref) |
| `Polynomials.domain`       |   | Should return a `Polynomials.Interval` instance|

As always, if the default implementation does not work or there are more efficient ways of implementing, feel free to overwrite functions from `common.jl` for your type.

The general idea is the container type should provide the vector operations of polynomial addition, subtraction, and scalar multiplication.
The latter is generically implemented through a `map(f,p)` method. The second example illustrates, though it isn't expected that container types will need being defined by users of this package.

The basis type directs dispatch for other operations and allows definitions for `one` and `variable`. An `evalpoly` method may be defined for a given basis type, though specializations based on the container may be desirable.

Methods like `*` will typically need to consider both the underlying container type and the basis, though if `convert` methods are defined, the defaults can be utilized as converting to the `Polynomial` type, performing the operation, then converting back is possible, though likely not as efficient.


!!! note
    Most promotion rules will coerce towards the [`Polynomial`](@ref) type, so not all methods have to be implemented if you provide a conversion function.



## A new basis type

The generalized Laguerre polynomials are orthogonal polynomials parameterized  by ``\alpha`` and defined recursively by

```math
\begin{align*}
L^\alpha_1(x) &= 1\\
L^\alpha_2(x) &= 1 + \alpha - x\\
L^\alpha_{n+1}(x) &= \frac{2n+1+\alpha -x}{n+1} L^\alpha_n(x) - \frac{n+\alpha}{n+1} L^\alpha_{n-1}(x)\\
&= (A_nx +B_n) \cdot L^\alpha_n(x) - C_n \cdot L^\alpha_{n-1}(x).
\end{align*}
```

There are other [characterizations available](https://en.wikipedia.org/wiki/Laguerre_polynomials). The three-point recursion, described by `A`,`B`, and `C` is used below for evaluation.

We show how to define a new basis type, `LaguerreBasis`, leveraging one of the existing container types.
In this example our basis type has a parameter. The `ChebyshevT` type, gives a related example of how this task can be implemented.


First we load the package and import a few non-exported functions:

```jldoctest abstract_univariate_polynomial
julia> using Polynomials;

julia> import Polynomials: AbstractUnivariatePolynomial, AbstractBasis, MutableDensePolynomial;
```

We define the basis with:

```jldoctest abstract_univariate_polynomial
julia> struct LaguerreBasis{alpha} <: AbstractBasis end

julia> Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{LaguerreBasis{α},T,X}}) where {α,T,X} =
           "L^$(α)"
```

We added a method to `basis_symbol` to show this basis. The display of the basis symbol has a poor default. The method above requires the full type, as the indeterminate, `X`, may be part of the desired output.  More generally, `Polynomials.printbasis` can have methods added to adjust for different display types.

Polynomial types can be initiated through specifying a storage type and a basis type, say:

```jldoctest abstract_univariate_polynomial
julia> P = MutableDensePolynomial{LaguerreBasis{0}}
MutableDensePolynomial{LaguerreBasis{0}}
```

Instances can now be created:

```jldoctest abstract_univariate_polynomial
julia> p = P([1,2,3])
MutableDensePolynomial(1L^0_0 + 2*L^0_1 + 3*L^0_2)
```

Or using other storage types:

```jldoctest abstract_univariate_polynomial
julia> Polynomials.ImmutableDensePolynomial{LaguerreBasis{1}}((1,2,3))
Polynomials.ImmutableDensePolynomial(1L^1_0 + 2*L^1_1 + 3*L^1_2)
```

All polynomial types have vector addition and scalar multiplication defined, as these are basis independent:

```jldoctest abstract_univariate_polynomial
julia> q = P([1,2])
MutableDensePolynomial(1L^0_0 + 2*L^0_1)

julia> p + q
MutableDensePolynomial(2L^0_0 + 4*L^0_1 + 3*L^0_2)

julia> 2p
MutableDensePolynomial(2L^0_0 + 4*L^0_1 + 6*L^0_2)
```

For a new basis, there are no default methods for polynomial evaluation and polynomial multiplication; and no defaults for `one` (used by default for scalar addition), and `variable` (used by default in conversion).

For the Laguerre Polynomials, Clenshaw recursion can be used for evaluation.

```jldoctest abstract_univariate_polynomial
julia> function ABC(::Type{LaguerreBasis{α}}, n) where {α}
           o = one(α)
           d = n + o
           (A=-o/d, B=(2n + o + α)/d, C=(n+α)/d)
       end
ABC (generic function with 1 method)
```

```jldoctest abstract_univariate_polynomial
julia> function clenshaw_eval(p::P, x::S) where {α, Bᵅ<: LaguerreBasis{α}, T, P<:AbstractUnivariatePolynomial{Bᵅ,T}, S}
           d = degree(p)
           R = typeof(((one(α) * one(T)) * one(S)) / 1)
           p₀ = one(R)
           d == -1 && return zero(R)
           d == 0 && return p[0] * one(R)
           Δ0 = p[d-1]
           Δ1 = p[d]
           @inbounds for i in (d - 1):-1:1
               A,B,C = ABC(Bᵅ, i)
               Δ0, Δ1 =
                   p[i] - Δ1 * C, Δ0 + Δ1 * muladd(x, A, B)
           end
           A,B,C = ABC(Bᵅ, 0)
           p₁ = muladd(x, A, B) * p₀
           return Δ0 * p₀ + Δ1 * p₁
       end
clenshaw_eval (generic function with 1 method)
```
Internally, `evalpoly` is called so we forward that method.

```jldoctest abstract_univariate_polynomial
julia> Polynomials.evalpoly(x, p::P) where {P<:AbstractUnivariatePolynomial{<:LaguerreBasis}} =
               clenshaw_eval(p, x)
```

We test this out by passing in the variable `x` in the standard basis:

```jldoctest abstract_univariate_polynomial
julia> p = P([0,0,1])
MutableDensePolynomial(L^0_2)

julia> x = variable(Polynomial)
Polynomial(1.0*x)

julia> p(x)
Polynomial(1.0 - 2.0*x + 0.5*x^2)
```

This shows evaluation works and also that conversion to the `Polynomial` type is available through polynomial evaluation. This is used by default by `convert`, so we immediately have other `convert` methods available:

```jldoctest abstract_univariate_polynomial
julia> convert(ChebyshevT, p)
ChebyshevT(1.25⋅T_0(x) - 2.0⋅T_1(x) + 0.25⋅T_2(x))
```

Or, using some extra annotations to have rational arithmetic used, we can compare to easily found representations in the standard basis:

```jldoctest abstract_univariate_polynomial
julia> q = Polynomials.basis(MutableDensePolynomial{LaguerreBasis{0//1}, Int}, 5)
MutableDensePolynomial(L^0//1_5)

julia> x = variable(Polynomial{Int})
Polynomial(x)

julia> q(x)
Polynomial(1//1 - 5//1*x + 5//1*x^2 - 5//3*x^3 + 5//24*x^4 - 1//120*x^5)
```

The values of `one` and `variable` are straightforward to define, as ``L_0=1`` and ``L_1=1 - x`` or ``x = L_0 - L_1``

```jldoctest abstract_univariate_polynomial
julia> Polynomials.one(::Type{P}) where {B<:LaguerreBasis,T,X,P<:AbstractUnivariatePolynomial{B,T,X}} =
           P([one(T)])

julia> Polynomials.variable(::Type{P}) where {B<:LaguerreBasis,T,X,P<:AbstractUnivariatePolynomial{B,T,X}} =
           P([one(T), -one(T)])
```

To see this is correct, we have:

```jldoctest abstract_univariate_polynomial
julia> variable(P)(x) == x
true
```

Scalar addition defaults to a call to `one(p)`, so this is now defined:

```jldoctest abstract_univariate_polynomial
julia> 2 + p
MutableDensePolynomial(2L^0_0 + L^0_2)
```

Often it is more performant to implement a specific method for `scalar_add`. Here we utilize the fact that ``L_0 = 1`` to manipulate the coefficients. Below we specialize to a container type:

```jldoctest abstract_univariate_polynomial
julia> function Polynomials.scalar_add(c::S, p::P) where {B<:LaguerreBasis,T,X,
                                                          P<:MutableDensePolynomial{B,T,X},S}
           R = promote_type(T,S)
           iszero(p) && return MutableDensePolynomial{B,R,X}(c)
           cs = convert(Vector{R}, copy(p.coeffs))
           cs[1] += c
           MutableDensePolynomial{B,R,X}(cs)
       end

julia> p + 3
MutableDensePolynomial(3L^0_0 + L^0_2)
```

Multiplication defaults to a code path where the two polynomials are promoted to a common type and then multiplied.
Here we implement polynomial multiplication through conversion to the polynomial type. The [direct formula](https://londmathsoc.onlinelibrary.wiley.com/doi/pdf/10.1112/jlms/s1-36.1.399) could be implemented, but that isn't so illustrative for this example. See the `SpecialPolynomials` package for an implementation.

```jldoctest abstract_univariate_polynomial
julia> function Base.:*(p::MutableDensePolynomial{B,T,X},
                        q::MutableDensePolynomial{B,S,X}) where {B<:LaguerreBasis, T,S,X}
           x = variable(Polynomial{T,X})
           p(x) * q(x)
       end
```

Were it defined, a `convert` method from `Polynomial` to the `LaguerreBasis` could be used to implement multiplication, as we have defined a `variable` method.

## A new container type

This example shows how to make a new container type, though this should be unnecessary, given the current variety, there may be gains to be had (e.g. an immutable, sparse type?)
In this case, we offer a minimal example where the polynomial type aliases the vector defining the coefficients is created.  For other bases, more methods may be necessary to define (again, refer to ChebyshevT for an example).


We have two constructor methods. The first is the typical code path. It makes a copy of the coefficients and then wraps those within the polynomial container type. For performance reasons, generically it is helpful to pass in a flag to indicate no copying or checking of the input is needed (`Val{false}`). This is used by some inherited methods when we specialize to the `StandardBasis` type. Generically, a container type *may* accept an offset, though this type won't; a `0`-based vector is implicit.

```jldoctest new_container_type
julia> using Polynomials

julia> struct AliasPolynomialType{B,T,X} <: Polynomials.AbstractDenseUnivariatePolynomial{B, T, X}
           coeffs::Vector{T}
           function AliasPolynomialType{B, T, X}(coeffs::AbstractVector{S}, o::Int=0) where {B, T, S, X}
               new{B,T,Symbol(X)}(convert(Vector{T}, copy(coeffs)))
           end
           function AliasPolynomialType{B, T, X}(::Val{false}, coeffs::AbstractVector{S}, o::Int=0) where {B, T, S, X}
               new{B,T,Symbol(X)}(convert(Vector{T}, coeffs))
           end
       end

julia> Polynomials.@poly_register AliasPolynomialType
```

The call to `@poly_register` adds many different means to construct polynomials of this type along with some other default methods.

A few methods need defining to get indexing to work:

```jldoctest new_container_type
julia> Base.firstindex(p::AliasPolynomialType) = 0

julia> Base.lastindex(p::AliasPolynomialType) = length(p.coeffs) - 1

```

```jldoctest new_container_type
julia> Polynomials.constructorof(::Type{<:AliasPolynomialType{B}}) where {B} = AliasPolynomialType{B}

```

We need to add in the vector-space operations:

```jldoctest new_container_type
julia> function Base.:+(p::AliasPolynomialType{B,T,X}, q::AliasPolynomialType{B,S,X}) where {B,S,T,X}
           R = promote_type(T,S)
           n = maximum(degree, (p,q))
           cs = [p[i] + q[i] for i in 0:n]
           AliasPolynomialType{B,R,X}(Val(false), cs)  # save a copy
       end

julia> function Base.:-(p::AliasPolynomialType{B,T,X}, q::AliasPolynomialType{B,S,X}) where {B,S,T,X}
           R = promote_type(T,S)
           n = maximum(degree, (p,q))
           cs = [p[i] - q[i] for i in 0:n]
           AliasPolynomialType{B,R,X}(Val(false), cs)
       end

julia> function Base.map(fn, p::P) where {B,T,X,P<:AliasPolynomialType{B,T,X}}
           cs = map(fn, p.coeffs)
           R = eltype(cs)
           AliasPolynomialType{B,R,X}(Val(false), cs)
       end

```


A type and a basis defines a polynomial type.
This example uses the  `StandardBasis` basis type and consequently inherits the methods mentioned above that otherwise would need implementing.

```jldoctest new_container_type
julia> AliasPolynomial = AliasPolynomialType{Polynomials.StandardBasis};

```

To see this new polynomial type in action, we have:

```jldoctest new_container_type
julia> xs = [1,2,3,4];

julia> p = AliasPolynomial(xs)
AliasPolynomialType(1 + 2*x + 3*x^2 + 4*x^3)

julia> q = AliasPolynomial(1.0, :y)
AliasPolynomialType(1.0)

julia> 2p - q
AliasPolynomialType(3.0 + 4.0*x + 6.0*x^2 + 8.0*x^3)

julia> (derivative ∘ integrate)(p) == p
true

julia> p(3)
142
```

The default for polynomial multiplication is to call `*` for two instances of the type with the same variable, and possibly different element types. For standard basis types, we can add this method:

```jldoctest new_container_type
julia> Base.:*(p::AliasPolynomialType{T,X}, q::AliasPolynomialType{S,X}) where {T,S,X} = Polynomials._standard_basis_multiplication(p,q)

julia> p * p
AliasPolynomialType(1 + 4*x + 10*x^2 + 20*x^3 + 25*x^4 + 24*x^5 + 16*x^6)
```


For the Polynomial type, the default on operations is to copy the array. For this type, it might seem reasonable -- to avoid allocations -- to update the coefficients in place for scalar addition and scalar multiplication.

Scalar addition, `p+c`, defaults to `p + c*one(p)`, or polynomial addition, which is not inplace without additional work. As such, we create a new method and an infix operator

```jldoctest new_container_type
julia> function scalar_add!(c::T, p::AliasPolynomial{T}) where {T}
           p.coeffs[1] += c
           p
       end;

julia> p::AliasPolynomial +ₛ c::Number = scalar_add!(c, p);

julia> c::Number +ₛ p::AliasPolynomial = scalar_add!(c, p);
```

The viewpoint that a polynomial represents a vector of coefficients leads to an expectation that vector operations should match when possible. Scalar multiplication is a vector operation, so it seems reasonable to override the broadcast machinery to implement an in place operation (e.g. `p .*= 2`). By default, the polynomial types are not broadcastable over their coefficients. We would need to make a change there and modify the `copyto!` function:


```jldoctest new_container_type
julia> Base.broadcastable(p::AliasPolynomial) = p.coeffs;

julia> Base.ndims(::Type{<:AliasPolynomial}) = 1

julia> Base.copyto!(p::AliasPolynomial, x) = (copyto!(p.coeffs, x); chop!(p));

julia> p
AliasPolynomialType(1 + 2*x + 3*x^2 + 4*x^3)

julia> p .*= 2
AliasPolynomialType(2 + 4*x + 6*x^2 + 8*x^3)

julia> p ./= 2
AliasPolynomialType(1 + 2*x + 3*x^2 + 4*x^3)
```

Trying to divide again would throw an error, as the result would not fit with the integer type of `p`.

Now `p` is treated as the vector `p.coeffs`, as regards broadcasting, so some things may be surprising, for example this expression returns a vector, not a polynomial:

```jldoctest new_container_type
julia> p .+ 2
4-element Vector{Int64}:
 3
 4
 5
 6
```

The unexported `Polynomials.PnPolynomial` polynomial type implements much of the above.

----

```@docs
Polynomials.AbstractUnivariatePolynomial
Polynomials.AbstractBasis
```
