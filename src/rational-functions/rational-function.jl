"""
    RationalFunction(p::AbstractPolynomial, q::AbstractPolynomial)
    p // q

Create a rational expression (`p//q`) from the two polynomials. 

Common factors are not cancelled by the constructor, as they are for
the base `Rational` type. The [`lowest_terms(pq)`](@ref) function attempts
that operation.

For purposes of iteration, a rational function is treated like a two-element container.

## Examples
```
julia> using Polynomials

julia> p,q = fromroots(Polynomial, [1,2,3]), fromroots(Polynomial, [2,3,4])
(Polynomial(-6 + 11*x - 6*x^2 + x^3), Polynomial(-24 + 26*x - 9*x^2 + x^3))

julia> pq = p // q
(-6 + 11*x - 6*x^2 + x^3) // (-24 + 26*x - 9*x^2 + x^3)

julia> lowest_terms(pq)
(-0.333333 + 0.333333*x) // (-1.33333 + 0.333333*x)

julia> pq(2.5)
-1.0

julia> pq(2) # uses first non-`0/0` ratio of `p⁽ᵏ⁾/q⁽ᵏ⁾`
-0.5

julia> pq^2
(36 - 132*x + 193*x^2 - 144*x^3 + 58*x^4 - 12*x^5 + x^6) // (576 - 1248*x + 1108*x^2 - 516*x^3 + 133*x^4 - 18*x^5 + x^6)

julia> derivative(pq)
(-108 + 180*x - 111*x^2 + 30*x^3 - 3*x^4) // (576 - 1248*x + 1108*x^2 - 516*x^3 + 133*x^4 - 18*x^5 + x^6)
```

!!! note
    The [RationalFunctions.jl](https://github.com/aytekinar/RationalFunctions.jl) package was a helpful source of ideas. 

!!! note
    The `ImmutablePolynomial` type can not be used for rational functions, as the type requires the numerator and denominator to have the exact same type.


"""
struct RationalFunction{T, X, P<:AbstractPolynomial{T,X}} <: AbstractRationalFunction{T,X,P}
    num::P
    den::P
    function RationalFunction(p::P, q::P) where {T,X, P<:AbstractPolynomial{T,X}}
        new{T,X,P}(p, q)
    end
    function RationalFunction(p::P, q::T) where {T,X, P<:AbstractPolynomial{T,X}}
        new{T,X,P}(p, q*one(P))
    end
    function RationalFunction(p::T, q::Q) where {T,X, Q<:AbstractPolynomial{T,X}}
        new{T,X,Q}(p*one(Q), q)
    end
end

RationalFunction(p,q)  = RationalFunction(promote(p,q)...)
RationalFunction(p::ImmutablePolynomial,q::ImmutablePolynomial) = throw(ArgumentError("Sorry, immutable #polynomials are not a valid polynomial type for RationalFunction"))
function RationalFunction(p::LaurentPolynomial,q::LaurentPolynomial)
    𝐩 = convert(RationalFunction, p)
    𝐪 = convert(RationalFunction, q)
    𝐩 // 𝐪
end
RationalFunction(p::LaurentPolynomial,q::Number) = convert(RationalFunction, p) // q
RationalFunction(p::Number,q::LaurentPolynomial) = q // convert(RationalFunction, p)

RationalFunction(p::AbstractPolynomial) = RationalFunction(p,one(p))

# evaluation
(pq::RationalFunction)(x) = eval_rationalfunction(x, pq)

# Look like rational numbers
function Base.://(p::AbstractPolynomial,q::AbstractPolynomial)
    RationalFunction(p,q)
end


