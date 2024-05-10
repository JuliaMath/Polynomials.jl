"""
    RationalFunction(p::AbstractPolynomial, q::AbstractPolynomial)
    p // q

Create a rational expression (`p//q`) from the two polynomials.

Common factors are not cancelled by the constructor, as they are for
the base `Rational` type. The [`lowest_terms`](@ref) function attempts
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
    function RationalFunction{T,X,P}(p, q) where {T,X, P<:AbstractPolynomial{T,X}}
        new{T,X,P}(p, q)
    end

end

function RationalFunction(p::P, q::P) where {T,X, P<:AbstractPolynomial{T,X}}
    RationalFunction{T,X,P}(p, q)
end

RationalFunction(p::AbstractPolynomial{T,X}, q::AbstractPolynomial{S,X}) where {T,S,X} =
    RationalFunction(promote(p,q)...)

function RationalFunction(p::P, q::T) where {T,X, P<:AbstractPolynomial{T,X}}
    RationalFunction(p, (q * one(p)))
end
function RationalFunction(p::T, q::Q) where {T,X, Q<:AbstractPolynomial{T,X}}
    RationalFunction(p * one(q),  q)
end

function RationalFunction(p::P,q::P) where {T, X, P <: LaurentPolynomial{T,X}}

    m,n = firstindex(p), firstindex(q)
    p′,q′ = _shift(p, -m), _shift(q, -n)
    if m-n ≥ 0
        return RationalFunction{T,X,P}(_shift(p′, m-n), q′)
    else
        return RationalFunction{T,X,P}(p′, _shift(q′, n-m))
    end
end

# RationalFunction(p,q)  = RationalFunction(convert(LaurentPolynomial,p), convert(LaurentPolynomial,q))


# special case Laurent
function lowest_terms(pq::PQ; method=:numerical, kwargs...) where {T,X,
                                                                   P<:LaurentPolynomial{T,X}, #StandardBasisPolynomial{T,X},
                                                                   PQ<:AbstractRationalFunction{T,X,P}}
    p,q = pqs(pq)
    p′,q′ = convert(Polynomial, p), convert(Polynomial,q)
    u,v,w = uvw(p′,q′; method=method, kwargs...)
    v′,w′ = convert(LaurentPolynomial, v), convert(LaurentPolynomial, w)
    rational_function(PQ, v′/w′[end], w′/w′[end])
end


RationalFunction(p::AbstractPolynomial) = RationalFunction(p,one(p))

# evaluation
(pq::RationalFunction)(x) = eval_rationalfunction(x, pq)

# Look like rational numbers
function Base.://(p::AbstractPolynomial,q::AbstractPolynomial)
    RationalFunction(p,q)
end


# promotion
function Base.promote(pq::RationalFunction{T,X,P},
                      rs::RationalFunction{S,X,Q}) where {T,S,X,P<:AbstractPolynomial{T,X},
                                                          Q<:AbstractPolynomial{S,X}}
    𝑃 = promote_type(P, Q)
    p,q = pq
    r,s = rs
    (convert(𝑃,p) // convert(𝑃,q), convert(𝑃,r) // convert(𝑃,s))
end
