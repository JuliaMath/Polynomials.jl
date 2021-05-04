export RationalFunction
export poles, residues
export lowest_terms

## ----
"""
    AbstractRationalFunction{T,X,P}

Abstract type for holding ratios of polynomials of type `P{T,X}`.

Default methods for basic arithmetic operations are provided.

Numeric methods to cancel common factors, compute the poles, and return the residues are provided.

!!! Note:
    Requires `VERSION >= v"1.2.0"`
"""
abstract type AbstractRationalFunction{T,X,P} end


function Base.show(io::IO, pq::AbstractRationalFunction)
    p,q = pqs(pq)
    print(io,"(")
    printpoly(io, p)
    print(io, ") // (")
    printpoly(io, q)
    print(io, ")")
end

## helper to make a rational function of given type
function rational_function(::Type{R}, p::P, q::Q) where {R<:AbstractRationalFunction,
                                                         T,X,   P<:AbstractPolynomial{T,X},
                                                         S,   Q<:AbstractPolynomial{S,X}}
    constructorof(R)(promote(p,q)...)
end


## ---- conversion
function Base.convert(::Type{PQ}, pq′::PQ′) where {T,X,P,PQ <: AbstractRationalFunction{T,X,P},
                                                   T′,X′,P′,PQ′<:AbstractRationalFunction{T′,X′,P′} }
    !isconstant(pq′) && assert_same_variable(X,X′)
    p′,q′=pqs(pq′)
    𝑷 = isconstant(pq′) ? P :  promote_type(P, P′)
    p,q = convert(𝑷, p′), convert(𝑷, q′)
    rational_function(PQ, p, q)
end


# R ⊂ R[x] ⊂ R(x)
function Base.convert(::Type{PQ}, p::Number) where {PQ <: AbstractRationalFunction}
   P = eltype(PQ)
   rational_function(PQ, p * one(P), one(P))
end

function Base.convert(::Type{PQ}, p::P) where {PQ <: AbstractRationalFunction, P<:AbstractPolynomial}
    Q = eltype(PQ)
    q = convert(Q, p)
    rational_function(PQ, q, one(q))
end

function Base.convert(::Type{P}, pq::PQ) where {P<:AbstractPolynomial, PQ<:AbstractRationalFunction}
    p,q = pqs(pq)
    isconstant(q) || throw(ArgumentError("Can't convert rational function with non-constant denominator to a polynomial."))
    convert(P, p) / constantterm(q)
end

function Base.convert(::Type{S}, pq::PQ) where {S<:Number, T,X,P,PQ<:AbstractRationalFunction}
    !isconstant(pq) && throw(ArgumentError("Can't convert non-constant rational function to a number"))
    S(pq(0))
end



# promotion rule to promote things upwards
function Base.promote_rule(::Type{PQ}, ::Type{PQ′}) where {T,X,P,PQ <: AbstractRationalFunction{T,X,P},
                                                           T′,X′,P′,PQ′<:AbstractRationalFunction{T′,X′,P′} }
    assert_same_variable(X,X′)
    PQ_, PQ′_ = constructorof(PQ), constructorof(PQ′)
    𝑷𝑸 = PQ_ == PQ′ ? PQ_ : RationalFunction
    𝑷 = constructorof(typeof(variable(P)+variable(P′)))
    𝑻 = promote_type(T,T′)
    𝑷𝑸{𝑻,X,𝑷{𝑻,X}}
end
Base.promote_rule(::Type{PQ}, ::Type{P}) where {PQ <: AbstractRationalFunction, P<:AbstractPolynomial} = PQ
function Base.promote_rule(::Type{PQ}, ::Type{S}) where {T,X, P<:AbstractPolynomial{T,X}, PQ <: AbstractRationalFunction{T,X,P}, S<:Number}
    R = promote_type(S,T)
    P′ = constructorof(P){R,X}
    constructorof(PQ){R,X,P′}
end


## Look like rational numbers
# The p//q constructor is reserved for the `RationalFunction` type, but these should all be defined
# for other possible types.
function Base.://(p::PQ,q::PQ′) where {PQ <: AbstractRationalFunction, PQ′ <: AbstractRationalFunction}
    p0,p1 = p
    q0,q1 = q
    rational_function(promote_type(PQ, PQ′), p0*q1, p1*q0)
end

function Base.://(p::Union{Number,AbstractPolynomial},q::PQ) where {PQ <: AbstractRationalFunction}
    q0,q1 = q
    rational_function(PQ, p*q1, q0)
end
function Base.://(p::PQ, q::Union{Number,AbstractPolynomial}) where {PQ <: AbstractRationalFunction}
    p0, p1 = p
    rational_function(PQ, p0, p1*q)
end

Base.://(p::AbstractPolynomial,q::Number) = p // (q*one(p))
Base.://(p::Number, q::AbstractPolynomial) = (p*one(q)) // q
    
    
function Base.copy(pq::PQ) where {PQ <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(PQ, p, q)
end


## ----

# requires these field names
Base.numerator(pq::AbstractRationalFunction) = pq.num
Base.denominator(pq::AbstractRationalFunction) = pq.den

# Treat a RationalFunction as a tuple (num=p, den=q)
Base.length(pq::AbstractRationalFunction) = 2
function Base.iterate(pq::AbstractRationalFunction, state=nothing)
    state == nothing && return (numerator(pq), 1)
    state == 1 && return (denominator(pq), 2)
    nothing
end
Base.collect(pq::AbstractRationalFunction{T,X,P}) where {T,X,P} = collect(P, pq)
Base.broadcastable(pq::AbstractRationalFunction) = Ref(pq)

Base.eltype(pq::Type{<:AbstractRationalFunction{T,X,P}}) where {T,X,P} = P
Base.eltype(pq::Type{<:AbstractRationalFunction{T,X}}) where {T,X} = Polynomial{T,X}
Base.eltype(pq::Type{<:AbstractRationalFunction{T}}) where {T} = Polynomial{T,:x}
Base.eltype(pq::Type{<:AbstractRationalFunction}) = Polynomial{Float64,:x}


"""
    pqs(pq) 

Return `(p,q)`, where `pq=p/q`, as polynomials.
"""
pqs(pq::AbstractRationalFunction) = (numerator(pq), denominator(pq))

## ----

Base.size(F::AbstractRationalFunction) = ()
function Base.isinf(pq::AbstractRationalFunction)
    p,q=pqs(pq)
    iszero(q) && !iszero(p)
end
function Base.isnan(pq::AbstractRationalFunction) 
    p,q= pqs(pq)
    iszero(p) && iszero(q)
end
function Base.iszero(pq::AbstractRationalFunction)
    p,q= pqs(pq)
    iszero(p) && !iszero(q)
end
function Base.isone(pq::AbstractRationalFunction)
    p,q = pqs(pq)
    isconstant(p) && isconstant(q) && p == q
end

                                                             

## -----

_indeterminate(::Type{<:AbstractRationalFunction}) = nothing
_indeterminate(::Type{PQ}) where {T,X,PQ<:AbstractRationalFunction{T,X}} = X
indeterminate(pq::PQ) where {PQ<:AbstractRationalFunction} = indeterminate(PQ)
function indeterminate(::Type{PQ}, var=:x) where {PQ<:AbstractRationalFunction}
    X′ = _indeterminate(PQ)
    X = X′ == nothing ? Symbol(var) : X′
    X
end

function isconstant(pq::AbstractRationalFunction; kwargs...)
    p,q = pqs(lowest_terms(pq, kwargs...))
    isconstant(p) && isconstant(q)
end
isconstant(::Number) = true

function constantterm(pq::AbstractRationalFunction; kwargs...)
    p,q = pqs(pq)
    isconstant(pq) && return constantterm(p)/constantterm(q)
    throw(ArgumentError("No constant term defined for a non-constant polynomial"))
end

function Base.zero(pq::R) where {R <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(R, zero(p), one(q))
end
function Base.one(pq::R) where {R <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(R, one(p), one(q))
end
function variable(pq::R) where {R <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(R, variable(p), one(q))
end

function variable(::Type{PQ}) where {PQ <: AbstractRationalFunction}
    x = variable(eltype(PQ))
    rational_function(PQ, x, one(x))
end

    
# use degree as largest degree of p,q after reduction
function degree(pq::AbstractRationalFunction)
    pq′ = lowest_terms(pq)
    maximum(degree.(pqs(pq′)))
end

# Evaluation
function eval_rationalfunction(x, pq::AbstractRationalFunction{T}) where {T}
    num, den = pqs(pq)
    dn, dd = degree(num), degree(den)
    md = min(dn, dd)
    result = num(x)/den(x)
    md < 0 && return result
    while md >= 0
        !isnan(result) && return result
        num,den = derivative(num), derivative(den)
        result = num(x)/den(x)
        md -= 1
    end

    x*NaN
end

# equality
import Base: ==
function ==(p::AbstractRationalFunction{T,X,P}, q::AbstractRationalFunction{S,Y,Q}) where {T,X,P,S,Y,Q}
    isconstant(p) && isconstant(q) && p(0) == q(0) && return true
    X == Y || return false
    p₀, p₁ = pqs(p)
    q₀, q₁ = pqs(q)
    p₀ * q₁ == q₀ * p₁ || return false
end

function ==(p::AbstractRationalFunction{T,X,P}, q::Union{AbstractPolynomial, Number}) where {T,X,P}
    ==(promote(p,q)...)
end

function ==(p::Union{AbstractPolynomial, Number}, q::AbstractRationalFunction{T,X,P}) where {T,X,P}
   ==(promote(p,q)...)
end


function Base.isapprox(pq₁::PQ₁, pq₂::PQ₂,
                  rtol::Real = sqrt(eps(float(real(promote_type(T,S))))),
                  atol::Real = zero(float(real(promote_type(T,S))))) where {T,X,P,PQ₁<:AbstractRationalFunction{T,X,P},
                                                                     S,Y,Q,PQ₂<:AbstractRationalFunction{S,Y,Q}}

    p₁,q₁ = pqs(pq₁)
    p₂,q₂ = pqs(pq₂)

    isapprox(p₁*q₂, q₁*p₂; rtol=rtol, atol=atol)
end


# Arithmetic
function Base.:-(pq::PQ) where {PQ <: AbstractRationalFunction}
    p, q = copy(pq)
    rational_function(PQ, -p, q)
end

Base.:+(p::Number, q::AbstractRationalFunction) = q + p
Base.:+(p::AbstractRationalFunction,  q::Number) = p + q*one(p)
Base.:+(p::AbstractPolynomial, q::AbstractRationalFunction) = q + p
Base.:+(p::AbstractRationalFunction,  q::AbstractPolynomial) = p + (q//one(q))
Base.:+(p::AbstractRationalFunction, q::AbstractRationalFunction) = sum(promote(p,q))
# type should implement this
function Base.:+(p::R, q::R) where {T,X,P,R <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(R, p0*q1 + p1*q0, p1*q1)
end

Base.:-(p::Number, q::AbstractRationalFunction) = -q +  p
Base.:-(p::AbstractRationalFunction,  q::Number) = p - q*one(p)
Base.:-(p::AbstractPolynomial, q::AbstractRationalFunction) = -q + p
Base.:-(p::PQ,  q::AbstractPolynomial) where {PQ <: AbstractRationalFunction} = p - rational_function(PQ,q, one(q))
function Base.:-(p::AbstractRationalFunction, q::AbstractRationalFunction)
    p′, q′ = promote(p,q)
    p′ - q′
end
# type should implement this
function Base.:-(p::R, q::R) where {T,X,P,R <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(R, p0*q1 - p1*q0, p1*q1)
end

function Base.:*(p::Number, q::R) where {T, X, R <: AbstractRationalFunction{T,X}}
    q0,q1 = pqs(q)
    rational_function(R, (p*q0), q1)
end
function Base.:*(p::R,  q::Number) where {R <:AbstractRationalFunction}
    p0,p1 = pqs(p)
    rational_function(R, p0*q, p1)
end
function Base.:*(p::AbstractPolynomial, q::R) where {R <: AbstractRationalFunction}
    rational_function(R, p, one(p)) * q
end
Base.:*(p::R,  q::AbstractPolynomial) where {R <: AbstractRationalFunction} = p * rational_function(R,q, one(q))
# type should implement this
Base.:*(p::AbstractRationalFunction, q::AbstractRationalFunction) = prod(promote(p,q))
function Base.:*(p::R, q::R) where {T,X,P,R <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(R, p0*q0, p1*q1)
end



function Base.:/(p::Number, q::R) where {R <: AbstractRationalFunction}
    q0,q1 = pqs(q)
    rational_function(R,p*q1, q0)
end
function Base.:/(p::R,  q::Number) where {R <: AbstractRationalFunction}
    p0,p1 = pqs(p)
    rational_function(R, p0, (p1*q))
end
Base.:/(p::AbstractPolynomial, q::PQ) where {PQ <: AbstractRationalFunction} = rational_function(PQ, p,one(p)) / q
function Base.:/(p::PQ,  q::AbstractPolynomial) where {PQ <: AbstractRationalFunction} 
    p0,p1 = pqs(p)
    rational_function(PQ,p0, p1*q)
end
function Base.:/(p::AbstractRationalFunction, q::AbstractRationalFunction)
    p′,q′ = promote(p,q)
    p′ / q′
end
# type should implement this
function Base.:/(p::PQ, q::PQ) where {T,X,P,PQ <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(PQ, p0*q1, p1*q0)
end

function Base.:^(pq::P, n::Int) where {P <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(P, p^n, q^n)
end

function Base.inv(pq::P) where {P <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(P, q, p)
end

# conj, transpose... TODO

## derivative and integrals
function derivative(pq::P, n::Int=1) where {P <: AbstractRationalFunction}
    n <= 0 && return pq
    while n >= 1
        p,q = pqs(pq)
        pq = rational_function(P, (derivative(p)*q - p * derivative(q)), q^2)
        n -= 1
    end
    pq
end

function integrate(pq::P) where {P <: AbstractRationalFunction}
    p,q = pqs(pq)
    isconstant(q) && return rational_function(P, integrate(p), q)
    # XXX could work here, e.g.:
    # d,r = partial_fraction
    # ∫d + Σ∫r for each residue (issue with logs)
    throw(ArgumentError("Can only integrate rational functions with constant denominators"))
end

## ----

"""
    divrem(pq::AbstractRationalFunction; method=:numerical, kargs...)

Return `d,r` with `p/q = d + r/q` where `degree(numerator(r)) < degree(denominator(q))`, `d` a Polynomial, `r` a `AbstractRationalFunction`.

* `method`: passed to `gcd`
* `kwargs...`: passed to `gcd`
"""
function Base.divrem(pq::PQ; method=:numerical, kwargs...) where {PQ <: AbstractRationalFunction}
    p,q = pqs(pq)
    degree(p) < degree(q) && return (zero(p), pq)

    d,r = divrem(p,q)
    (d, rational_function(PQ, r, q))
    
end


# like Base.divgcd in rational.jl
# divide p,q by u
function _divgcd(V::Val{:euclidean}, pq; kwargs...)
    p, q = pqs(pq)
    u = gcd(V,p,q; kwargs...)
    p÷u, q÷u
end
function _divgcd(V::Val{:noda_sasaki}, pq; kwargs...)
    p, q = pqs(pq)
    u = gcd(V,p,q; kwargs...)
    p÷u, q÷u
end
function _divgcd(v::Val{:numerical}, pq; kwargs...)
    u,v,w,θ,κ = ngcd(pqs(pq)...; kwargs...) # u⋅v=p, u⋅w=q
    v, w
end


"""
    lowest_terms(pq::AbstractRationalFunction, method=:numerical)
    
Find GCD of `(p,q)`, `u`, and return `(p÷u)//(q÷u)`. Commonly referred to as lowest terms.
    
* `method`: passed to `gcd(p,q)`
* `kwargs`: passed to `gcd(p,q)`

By default, `AbstractRationalFunction` types do not cancel common factors. This method will numerically cancel common factors, returning the normal form, canonicalized here by `q[end]=1`. The result and original may be considered equivalent as rational expressions, but different when seen as functions of the indeterminate.

"""
function lowest_terms(pq::PQ; method=:numerical, kwargs...) where {T,X,
                                                                   P<:StandardBasisPolynomial{T,X},
                                                                   PQ<:AbstractRationalFunction{T,X,P}}
    v,w = _divgcd(Val(method), pq; kwargs...)
    rational_function(PQ, v/w[end], w/w[end])
end

## ---- zeros, poles, ...
"""
    poles(pq::AbstractRationalFunction; method=:numerical, kwargs...)

For a rational function `p/q`, first reduces to normal form, then finds the roots and multiplicities of the resulting denominator.

"""
function poles(pq::AbstractRationalFunction; method=:numerical,  kwargs...)
    pq′ = lowest_terms(pq; method=method, kwargs...)
    den = denominator(pq′)
    mr = Multroot.multroot(den)
    (zs=mr.values, multiplicities = mr.multiplicities)
end

"""
    roots(pq::AbstractRationalFunction; kwargs...)

Return the `zeros` of the rational function (after cancelling commong factors, the `zeros` are the roots of the numerator.

"""
function roots(pq::AbstractRationalFunction; method=:numerical,  kwargs...)
    pq′ = lowest_terms(pq; method=method, kwargs...)
    den = numerator(pq′)
    mr = Multroot.multroot(den)
    (zs=mr.values, multiplicities = mr.multiplicities)
end

"""
    residues(pq::AbstractRationalFunction; method=:numerical,  kwargs...)

If `p/q =d + r/q`, returns `d` and the residues of a rational fraction `r/q`.

First expresses `p/q =d + r/q` with `r` of lower degree than `q` through `divrem`. 
Then finds the poles of `r/q`.
For a pole, `λj` of multiplicity `k` there are `k` residues, 
`rⱼ[k]/(z-λⱼ)^k`, `rⱼ[k-1]/(z-λⱼ)^(k-1)`, `rⱼ[k-2]/(z-λⱼ)^(k-2)`, …, `rⱼ[1]/(z-λⱼ)`.
The residues are found using this formula: 
`1/j! * dʲ/dsʲ (F(s)(s - λⱼ)^k` evaluated at `λⱼ` ([5-28](https://stanford.edu/~boyd/ee102/rational.pdf)).


## Example

(From page 5-33 of above pdf)

```jldoctest
julia> s = variable(Polynomial, :s)
Polynomial(1.0*s)

julia> pq = (-s^2 + s + 1) // (s * (s+1)^2)
(1.0 + 1.0*s - 1.0*s^2) // (1.0*s + 2.0*s^2 + 1.0*s^3)

julia> d,r = residues(pq)
(Polynomial(0.0), Dict(-1.0 => [-1.9999999999999978, 1.0000000000000029], 2.0153919257182735e-18 => [1.0]))

julia> iszero(d)
true

julia> z = variable(pq)
(1.0*s) // (1.0)

julia> for (λ, rs) ∈ r # reconstruct p/q from output of `residues`
           for (i,rᵢ) ∈ enumerate(rs)
               d += rᵢ/(z-λ)^i
           end
       end

julia> p′, q′ = lowest_terms(d)
(1.0 + 1.0*s - 1.0*s^2) // (6.64475e-16 + 1.0*s + 2.0*s^2 + 1.0*s^3)

julia> q′ ≈ (s * (s+1)^2) # works, as q is monic
true

julia> p′ ≈ (-s^2 + s + 1) 
true
```



!!! Note:
    There are several areas where numerical issues can arise. The `divrem`, the identification of multiple roots (`multroot`), the evaluation of the derivatives, ...

"""
function residues(pq::AbstractRationalFunction; method=:numerical,  kwargs...)

    
    d,r′ = divrem(pq)
    r = lowest_terms(r′; method=method, kwargs...)
    b,a = pqs(r)
    a′ = derivative(a)

    residues = Any[]
    mr = Multroot.multroot(a)

    for (λₖ, mₖ) ∈ zip(mr.values, mr.multiplicities)

        if mₖ == 1
            push!(residues, λₖ => [b(λₖ)/a′(λₖ)])
        else
            # returns rₖ₁, rₖ₂,...rₖₘ where rₖ,i/(s-λₖ)ⁱ is part of the decomposition
            s = variable(a)
            F = lowest_terms(r*(s-λₖ)^mₖ)
            rs = [F(λₖ)]
            j! = 1
            for j ∈ 1:mₖ-1
                dF = lowest_terms(derivative(F))
                pushfirst!(rs, 1/j! * dF(λₖ))
                j! *= (j+1)
            end
            push!(residues, λₖ => rs)
        end
    end
    d, Dict(residues...)
end

"""
    partial_fraction(pq::AbstractRationalFunction; method=:numerical, kwargs...)

For a rational function `p/q = d + r/q`, with `degree(r) < degree(q)` returns `d` and
the terms that comprise `r/q`: For each pole with multiplicity returns
`rⱼ[k]/(z-λⱼ)^k`, `rⱼ[k-1]/(z-λⱼ)^(k-1)`, `rⱼ[k-2]/(z-λⱼ)^(k-2)`, …, `rⱼ[1]/(z-λⱼ)`.

Should be if `p/q` is in normal form and `d,r=partial_fraction(p//q)` that
`d + sum(r) - p//q ≈ 0`

"""
function partial_fraction(pq::AbstractRationalFunction; method=:numerical, kwargs...)
    d,r = residues(pq; method=method, kwargs...)
    s = variable(pq)
    d, partial_fraction(Val(:residue), r, s)
end

function partial_fraction(::Val{:residue}, r, s::PQ) where {PQ}
    terms = []
    for (λₖ,rsₖ) ∈ r
        for  (rⱼ,j) ∈ zip(rsₖ, length(rsₖ):-1:1)
            push!(terms, rⱼ/(s-λₖ)^j)
        end
    end
    terms

end
            
            
            
            
    

    





