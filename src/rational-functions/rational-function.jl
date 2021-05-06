"""
    RationalFunction(p::AbstractPolynomial, q::AbstractPolynomial)
    p // q

Create a rational expression (`p/q`) from the two polynomials. 

There is no attempt to cancel common factors. The [`lowest_terms`](@ref) function attempts to do that.

For purposes of iteration, a rational function is treated like a tuple.

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

!!! Note:
    The [RationalFunctions.jl](https://github.com/aytekinar/RationalFunctions.jl) was a helpful source of ideas.


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
RationalFunction(p::AbstractPolynomial) = RationalFunction(p,one(p))

# evaluation
(pq::RationalFunction)(x) = eval_rationalfunction(x, pq)

# Look like rational numbers
function Base.://(p::AbstractPolynomial,q::AbstractPolynomial)
    RationalFunction(p,q)
end


