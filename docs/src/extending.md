# Extending Polynomials

The [`AbstractUnivaeriatePolynomial`](@ref) type was made to be extended.

```@docs
Polynomials.AbstractUnivariatePolynomial
```

A polynomial's  coefficients  are  relative to some *basis*. The `Polynomial` type relates coefficients  `[a0, a1,  ..., an]`, say,  to the  polynomial  `a0 +  a1*x + a2*x^  + ... +  an*x^n`,  through the standard  basis  `1,  x,  x^2, ..., x^n`.  New polynomial  types typically represent the polynomial through a different  basis. For example,  `CheyshevT` uses a basis  `T_0=1, T_1=x,  T_2=2x^2-1,  ...,  T_n  =  2xT_{n-1} - T_{n-2}`.  For this type  the  coefficients  `[a0,a1,...,an]` are associated with  the polynomial  `a0*T0  + a1*T_1 +  ...  +  an*T_n`.
An `AbstractUnivariatePolynomial` polynomial consists of a basis and a storage type. There a a few different storage types.

A basis inherits from `Polynomials.AbstractBasis`, in the example our basis type has a parameter. The `ChebyshevT` type, gives a related example of how this task can be implemented.

To implement a new polynomial type, `P`, the following methods should be implemented:


| Function | Required | Notes |
|----------|:--------:|:------------|
| Constructor | x | |
| Type function (`(::P)(x)`) | x | |
| `one`| | Convenience to find constant in new basis |
| `variable`| | Convenience to find monomial `x` in new  basis|
| `Base.evalpoly(x, p::P)` | x |to evaluate the polynomial at `x` |
| `scalar_add(c::S, ::P)` | x | Scalar addition |
| `*(::P, ::P)` | | Multiplication of polynomials |
| `divrem` | | Required for [`gcd`](@ref)|
| `Polynomials.domain` | x | Should return a `Polynomials.Interval` instance|
| `vander` | | Required for [`fit`](@ref) |
| `companion` | | Required for [`roots`](@ref) |
| `convert(::Polynomial, ...)` | | Defaults to polynomial evaluation |
| `convert(::P, p::Polynomial)`| | Default of polynomial evaluation requires `*` to be defined.|

Check out [`ChebyshevT`](@ref) for examples of this interface being extended.

As always, if the default implementation does not work or there are more efficient ways of implementing, feel free to overwrite functions from `common.jl` for your type.

!!! note
    Promotion rules will coerce towards the [`Polynomial`](@ref) type, so not all methods have to be implemented if you provide a conversion function.



### The generalized Laguerre polynomials

These are orthogonal polynomials parameterized  by $\alpha$ and defined recursively by

```math
\begin{align*}
L^\alpha_1(x) &= 1\\
L^\alpha_2(x) &= 1 + \alpha - x\\
L^\alpha_{n+1}(x) &= \frac{2n+1+\alpha -x}{n+1} L^\alpha_n(x) - \frac{n+\alpha}{n+1} L^\alpha_{n-1}(x)\\
&= (A_nx +B_n) \cdot L^\alpha_n(x) - C_n \cdot L^\alpha_{n-1}(x).
\end{align*}
```

There are other [characterizations available](https://en.wikipedia.org/wiki/Laguerre_polynomials). The three-point recursion, described by `A`,`B`, and `C` is used below for evaluation.

We define the basis with

```jldoctest abstract_univariate_polynomial
julia> using Polynomials;

julia> import Polynomials: AbstractUnivariatePolynomial, AbstractBasis, MutableDensePolynomial;

julia> struct LaguerreBasis{alpha} <: AbstractBasis end

julia> Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{LaguerreBasis{α},T,X}}) where {α,T,X} =
           "L^$(α)"
```

The basis symbol has a poor default. The method requires the full type, as the indeterminate, `X`, may part of the desired output. We added a method to `basis_symbol` to show this basis. More generally, `Polynomials.printbasis` can have methods added to adjust for different display types.

Polynomials can be initiated through specifying a storage type and a basis, say:

```jldoctest abstract_univariate_polynomial
julia> P = MutableDensePolynomial{LaguerreBasis{0}}
MutableDensePolynomial{LaguerreBasis{0}}

julia> p = P([1,2,3])
MutableDensePolynomial(1L^0_0 + 2*L^0_1 + 3*L^0_2)
```

Or using other storage types:

```jldoctest abstract_univariate_polynomial
julia> Polynomials.ImmutableDensePolynomial{LaguerreBasis{1}}((1,2,3))
Polynomials.ImmutableDensePolynomial(1L^1_0 + 2*L^1_1 + 3*L^1_2)
```

All polynomials have vector addition and scalar multiplication defined:

```jldoctest abstract_univariate_polynomial
julia> q = P([1,2])
MutableDensePolynomial(1L^0_0 + 2*L^0_1)

julia> p + q
MutableDensePolynomial(2L^0_0 + 4*L^0_1 + 3*L^0_2)
```

```jldoctest abstract_univariate_polynomial
julia> 2p
MutableDensePolynomial(2L^0_0 + 4*L^0_1 + 6*L^0_2)
```

For a new basis, there are no default methods for polynomial evaluation, scalar addition, and polynomial multiplication; and no defaults for `one`, and `variable`.

For the Laguerre Polynomials, Clenshaw recursion can be used for evaluation. Internally, `evalpoly` is called so we forward that method.

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

```jldoctest abstract_univariate_polynomial
julia> Polynomials.evalpoly(x, p::P) where {P<:AbstractUnivariatePolynomial{<:LaguerreBasis}} =
               clenshaw_eval(p, x)
```

We test it out by passing in the variable `x` in the standard basis:

```jldoctest abstract_univariate_polynomial
julia> p = P([0,0,1])
MutableDensePolynomial(L^0_2)

julia> x = variable(Polynomial)
Polynomial(1.0*x)

julia> p(x)
Polynomial(1.0 - 2.0*x + 0.5*x^2)
```

We see that conversion to the `Polynomial` type is available through polynomial evaluation. This is used by default, so we have `convert` methods available:

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


To implement scalar addition, we utilize the fact that ``L_0 = 1`` to manipulate the coefficients. Below we specialize to a container type:

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

The values of `one` and `variable` are straightforward, as ``L_0=1`` and ``L_1=1 - x`` or ``x = 1 - L_1``

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

Finally, we implement polynomial multiplication through conversion to the polynomial type. The [direct formula](https://londmathsoc.onlinelibrary.wiley.com/doi/pdf/10.1112/jlms/s1-36.1.399) could be implemented.

```jldoctest abstract_univariate_polynomial
julia> function Base.:*(p::MutableDensePolynomial{B,T,X},
                        q::MutableDensePolynomial{B,S,X}) where {B<:LaguerreBasis, T,S,X}
           x = variable(Polynomial{T,X})
           p(x) * q(x)
       end
```

Were it defined, a `convert` method from `Polynomial` to the `LaguerreBasis` could be used to implement multiplication.
