module RationalFunctionFit
using ..Polynomials
import ..Polynomials: RationalFunction, indeterminate, constructorof
import ..Polynomials: evalpoly
using LinearAlgebra

"""
    fit(::Type{RationalFunction}, xs::AbstractVector{S}, ys::AbstractVector{T}, m, n; var=:x)

Fit a rational function of the form `pq = (a₀ + a₁x¹ + … + aₘxᵐ) / (1 + b₁x¹ + … + bₙxⁿ)` to the data `(x,y)`.



!!! note
    This uses a simple implementation of the Gauss-Newton method
    to solve the non-linear least squares problem:
    `minᵦ Σ(yᵢ - pq(xᵢ,β)²`, where `β=(a₀,a₁,…,aₘ,b₁,…,bₙ)`.

    A more rapidly convergent method is used in the `LsqFit.jl`
    package, and if performance is important, re-expressing the
    problem for use with that package is suggested.

    Further, if an accurate rational function fit of adaptive degrees
    is of interest, the `BaryRational.jl` package provides an
    implementation of the AAA algorithm ("which offers speed,
    flexibility, and robustness we have not seen in other algorithms"
    [Nakatsukasa, Sète,
    Trefethen](https://arxiv.org/pdf/1612.00337.pdf)) and one using
    Floater-Hormann weights [Floater,
    Hormann](https://doi.org/10.1007/s00211-007-0093-y) ("that have no
    real poles and arbitrarily high approximation orders on any real
    interval, regardless of the distribution of the points")

    The [RationalApproximations](https://github.com/billmclean/RationalApproximations) package also has implementations of the AAA algorithm.

    A python library, [polyrat](https://github.com/jeffrey-hokanson/polyrat), has implementations of other algorithms.

## Example
```
julia> x = variable(Polynomial{Float64})
Polynomial(1.0*x)

julia> pq = (1+x)//(1-x)
(1.0 + 1.0*x) // (1.0 - 1.0*x)

julia> xs = 2.0:.1:3;

julia> ys = pq.(xs);

julia> v = fit(RationalFunction, xs, ys, 2, 2)
(1.0 + 1.0*x - 6.82121e-13*x^2) // (1.0 - 1.0*x + 2.84217e-13*x^2)

julia> maximum(abs, v(x)-pq(x) for x ∈ 2.1:0.1:3.0)
1.06314956838105e-12

julia> using BaryRational

julia> u = aaa(xs,ys)
(::BaryRational.AAAapprox{Vector{Float64}}) (generic function with 1 method)

julia> maximum(abs, u(x)-pq(x) for x ∈ 2.1:0.1:3.0)
4.440892098500626e-16

julia> u(variable(pq)) # to see which polynomial is used
(2.68328 + 0.447214*x - 1.78885*x^2 + 0.447214*x^3) // (2.68328 - 4.91935*x + 2.68328*x^2 - 0.447214*x^3)
```

"""
function Polynomials.fit(::Type{PQ}, xs::AbstractVector{S}, ys::AbstractVector{T}, m, n; var=:x) where {T,S, PQ<:RationalFunction}

    β₁,β₂ = gauss_newton(collect(xs), convert(Vector{float(T)}, ys), m, n)
    P = eltype(PQ)
    T′ = isnothing(Polynomials._eltype(P)) ? eltype(β₁) : eltype(P)
    X = indeterminate(PQ, var)
    P′ = constructorof(P){T′,X}
    p = P′(β₁)
    q = P′(vcat(1, β₂))

    p // q
end


"""
    fit(::Type{RationalFunction}, r::Polynomial, m, n; var=:x)

Fit a Pade approximant ([`pade_fit`](@ref)) to `r`.

Examples:

```jldoctext
julia> using Polynomials, PolynomialRatios

julia> x = variable()
Polynomial(x)

julia> ex = 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 # Taylor polynomial for e^x
Polynomial(1.0 + 1.0*x + 0.5*x^2 + 0.16666666666666666*x^3 + 0.041666666666666664*x^4 + 0.008333333333333333*x^5)

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 1,1)(x) for x ∈ 0:.05:0.5)
0.017945395966538547

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 1,2)(x) for x ∈ 0:.05:0.5)
0.0016624471707165078

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 2,1)(x) for x ∈ 0:.05:0.5)
0.001278729299871717

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 2,2)(x) for x ∈ 0:.05:0.5)
7.262205147950951e-5
```
"""
function Polynomials.fit(::Type{RationalFunction},r::Polynomial, m::Integer, n::Integer;var=:x)
    p,q = pade_fit(r, m,n, var=var)
    p // q
end



## ---- Pade
## https://mathworld.wolfram.com/PadeApproximant.html
"""
    pade_fit(r::Polynomial, m,n)

For a polynomial `r` of degree `d ≥ m + n`, find a rational function `p/q` with
`degree(p) ≤ m`, `degree(q) ≤ n` and `q*r - p = x^{m+n+1}*s(x)` for some polynomial `s`.

This implementation sets up a system of equations to identify `p` and `q`.
"""
function pade_fit(p::Polynomial{T}, m::Integer, n::Integer; var=:x) where {T}
    d = degree(p)
    @assert (0 <= m) && (1 <= n) && (m + n <= d)

    # could be much more performant
    c = convert(LaurentPolynomial, p) # for better indexing
    cs = [c[m+j-i] for j ∈ 1:n, i ∈ 0:n]

    qs′ = cs[:, 2:end] \ cs[:,1]
    qs = vcat(1, -qs′)

    cs = [c[0 + j - i] for j ∈ 0:m, i∈0:n]
    ps = cs * qs

    Polynomial(ps, var), Polynomial(qs,var)
end


## ---- Least Squares
## avoiding dependency on another package, LsqFit

using LinearAlgebra

# return pair ((a₀,...,aₙ), (b₁,...,bₘ))
function initial_guess(xs::Vector{T}, ys::Vector{S}, n, m) where {T, S}
    # p(x) = a₀ + ... + aₙx^n
    # q(x) = 1 + b₁x + ... + bₘx^m = 1 + r(x)
    # yᵢ + yᵢ * r(xᵢ) = p(xᵢ)
    # yᵢ = p(xᵢ) - r(xᵢ) * yᵢ
    k = n+1+m
    A = zeros(T, k, k)
    A[:,1] .= one(T)
    xs′ = xs[1:k]
    xs′′ = copy(xs′)
    for i ∈ 1:n
        A[:,1+i] = xs′
        xs′ .*= xs′′
    end
    xs′ = -copy(xs′′) .* ys[1:k]
    for i ∈ 1:m
        A[:, 1+n+i] = xs′
        xs′ .*= xs′′
    end
    β = pinv(A) *  ys[1:k]

end

function make_model(n)
    (x, β) -> begin
        β₁, β₂ = β[1:n+1], β[n+2:end]
        evalpoly.(x,(β₁,)) ./ evalpoly.(x, (vcat(1, β₂),))
    end
end

function J!(Jᵣ, xs::Vector{T}, β, n) where {T}
    β₁, β₂ = β[1:n+1],β[n+2:end]
    ps = evalpoly.(xs, (β₁,))
    qs = evalpoly.(xs, (vcat(1, β₂),))
    λ = one(T) ./ qs
    for i ∈ eachindex(β₁)
        Jᵣ[:,1] = λ
        λ .*= xs
    end
    λ = xs .* ps ./ (qs .* qs)
    for i ∈ eachindex(β₂)
        Jᵣ[:, n+1+i]=λ
        λ .*= xs
    end
    nothing
end


function gauss_newton(xs, ys::Vector{T}, n, m, tol=sqrt(eps(T))) where {T}

    β = initial_guess(xs, ys, n, m)
    model = make_model(n)

    Jᵣ = zeros(T, length(xs), 1 + n + m)

    Δ = norm(ys, Inf) * tol

    ϵₘ = norm(ys - model(xs, β), Inf)
    βₘ = copy(β)

    no_steps = 0

    while no_steps < 25
        no_steps += 1

        r = ys - model(xs, β)
        ϵ = norm(r, Inf)
        ϵ < Δ && return (β[1:n+1], β[n+2:end])
        if ϵ < ϵₘ
            ϵₘ = ϵ
            βₘ .= β
        end
        J!(Jᵣ, xs, β, n)
        Δᵦ = pinv(Jᵣ' * Jᵣ) * (Jᵣ' * r)
        β .-= Δᵦ

    end

    @warn "no convergence; returning best fit of many steps"
    return (βₘ[1:n+1], βₘ[n+2:end])
end


end
