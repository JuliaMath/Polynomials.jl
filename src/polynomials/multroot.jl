module Multroot

export multroot

using ..Polynomials
using LinearAlgebra

"""
    multroot(p; verbose=false, kwargs...)

Use `multroot` algorithm of
[Zeng](https://www.ams.org/journals/mcom/2005-74-250/S0025-5718-04-01692-8/S0025-5718-04-01692-8.pdf)
to identify roots of polynomials with suspected multiplicities over
`Float64` values, typically.

Example:

```jldoctest
julia> using Polynomials

julia> p = fromroots([sqrt(2), sqrt(2), sqrt(2), 1, 1])
Polynomial(-2.8284271247461907 + 11.656854249492383*x - 19.07106781186548*x^2 + 15.485281374238573*x^3 - 6.242640687119286*x^4 + 1.0*x^5)

julia> roots(p)
5-element Vector{ComplexF64}:
  0.999999677417768 + 0.0im
 1.0000003225831504 + 0.0im
 1.4141705716005981 + 0.0im
 1.4142350577588885 - 3.72273772278647e-5im
 1.4142350577588885 + 3.72273772278647e-5im

julia> m = Polynomials.Multroot.multroot(p);

julia> Dict(m.values .=> m.multiplicities)
Dict{Float64, Int64} with 2 entries:
  1.41421 => 3
  1.0     => 2
```

The algorithm has two stages. First it uses `pejorative_manifold` to
identify the number of distinct roots and their multiplicities. This
uses the fact if `p=Π(x-zᵢ)ˡⁱ`, `u=gcd(p, p′)`, and `u⋅v=p` then
`v=Π(x-zi)` is square free and contains the roots of `p` and `u` holds
the multiplicity details, which are identified by recursive usage of
`ncgd`, which identifies `u` and `v` above even if numeric
uncertainties are present.

Second it uses `pejorative_root` to improve a set of initial guesses
for the roots under the assumption the multiplicity structure is
correct using a Newton iteration scheme.

The following tolerances, passed through to `pejorative_manifold` by
`kwargs...`, are all used in the first stage, to identify the
multiplicity structure:

* `θ`: the singular value threshold, set to `1e-8`. This is used by
  `Polynomials.ngcd` to assess if a matrix is rank deficient by
  comparing the smallest singular value to `θ ⋅ ‖p‖₂`.

* `ρ`: the initial residual tolerance, set to `1e-10`. This is passed
  to `Polynomials.ngcd`, the GCD finding algorithm as a measure for
  the absolute tolerance multiplied by `‖p‖₂`.

* `ϕ`: A scale factor, set to `100`. As the `ngcd` algorithm is called
  recursively, this allows the residual tolerance to scale up to match
  increasing numeric errors.

Returns a named tuple with

* `values`: the identified roots
* `multiplicities`: the corresponding multiplicities
* `κ`: the estimated condition number
* `ϵ`: the backward error, `‖p̃ - p‖_w`.

If the wrong multiplicity structure is identified in step 1, then
typically either `κ` or `ϵ` will be large. The estimated forward
error, `‖z̃ - z‖₂`, is bounded up to higher order terms by `κ ⋅ ϵ`.
This will be displayed if `verbose=true` is specified.

For polynomials of degree 20 or higher, it is often the case the `l`
is misidentified.

"""
function multroot(p::Polynomials.StandardBasisPolynomial{T}; verbose=false,
                  kwargs...) where {T}

    # degenerate case, constant
    degree(p) == 0 && return (values=T[], multiplicities=Int[], κ=NaN, ϵ=NaN)

    # degenerate case, all zeros
    if (nz = findfirst(!iszero, coeffs(p))) == length(coeffs(p))
        return (values=zeros(T,1), multiplicities=[nz-1], κ=NaN, ϵ=NaN)
    end

    z, l = pejorative_manifold(p; kwargs...)
    z̃ = pejorative_root(p, z, l)
    κ, ϵ = stats(p, z̃, l)

    verbose && show_stats(κ, ϵ)
    (values = z̃, multiplicities = l, κ = κ, ϵ = ϵ)

end

# The multiplicity structure, `l`, gives rise to a pejorative manifold `Πₗ = {Gₗ(z) | z∈ Cᵐ}`.
function pejorative_manifold(p::Polynomials.StandardBasisPolynomial{T};
                             θ = 1e-8,  # zero singular-value threshold
                             ρ = 1e-10, # initial residual tolerance
                             ϕ = 100    # residual growth factor
                             )  where {T}
    p₂ = norm(p, 2)
    u, v, w, ρⱼ, κ = Polynomials.ngcd(p, derivative(p);
                                      satol = θ * p₂, srtol = zero(real(T)),
                                      atol = ρ * p₂, rtol = zero(real(T))
    )

    zs = roots(v)
    nrts = length(zs)
    ls = ones(Int, nrts)
    λ = ϕ
    while !Polynomials.isconstant(u)

        nu₂ = norm(u,2)
        θ′ = θ * nu₂
        ρ′ = max(ρ, ϕ * ρⱼ)  # paper uses just latter
        u, v, w, ρⱼ, κ = Polynomials.ngcd(u, derivative(u),
                                          satol = θ′,
                                          atol = ρ′,
                                          minⱼ = degree(u) - nrts
                                          )

        rts = roots(v)
        nrts = length(rts)

        for z in rts
            tmp, ind = findmin(abs.(zs .- z))
            ls[ind] += 1
        end

    end

    zs, ls # estimate for roots, multiplicities

end

"""
    pejorative_root(p, zs, ls; kwargs...)

Find a *pejorative* *root* for `p` given multiplicity structure `ls` and initial guess `zs`.

The pejorative manifold for a multplicity structure `l` is denoted `{Gₗ(z) | z ∈ Cᵐ}`. A pejorative
root is a least squares minimizer of `F(z) = W ⋅ [Gₗ(z) - a]`. Here `a ~ (p_{n-1}, p_{n-2}, …, p_1, p_0) / p_n` and `W` are weights given by `min(1, 1/|aᵢ|)`. When `l` is the mathematically correct structure, then `F` will be `0` for a pejorative root. If `l` is not correct, then the backward error `‖p̃ - p‖_w` is typically large, where `p̃ = Π (x-z̃ᵢ)ˡⁱ`.

This follows Algorithm 1 of [Zeng](https://www.ams.org/journals/mcom/2005-74-250/S0025-5718-04-01692-8/S0025-5718-04-01692-8.pdf)
"""
function pejorative_root(p::Polynomials.StandardBasisPolynomial,
                         zs::Vector{S}, ls; kwargs...) where {S}
    ps = (reverse ∘ coeffs)(p)
    pejorative_root(ps, zs, ls; kwargs...)
end

# p is [1, a_{n-1}, a_{n-2}, ..., a_1, a_0]
function pejorative_root(p, zs::Vector{S}, ls;
                 τ = eps(float(real(S)))^(2/3),
                 maxsteps = 100, # linear or quadratic
                 verbose=false
                 ) where {S}

    ## Solve WJ Δz = W(Gl(z) - a)
    ## using weights min(1/|aᵢ|), i ≠ 1

    m,n = sum(ls), length(zs)

    # storage
    a = p[2:end]./p[1]     # a ~ (p[n-1], p[n-2], ..., p[0])/p[n]
    W = Diagonal([min(1, 1/abs(aᵢ)) for aᵢ in a])
    J = zeros(S, m, n)
    G = zeros(S, 1 + m)
    Δₖ = zeros(S, n)
    zₖs = copy(zs)  # updated

    cvg = false
    δₖ₀ = -Inf
    for ctr in 1:maxsteps

        evalJ!(J, zₖs, ls)
        evalG!(G, zₖs, ls)
        Δₖ .= (W*J) \ (W*(view(G, 2:1+m) .- a)) # weighted least squares

        δₖ₁ = norm(Δₖ, 2)
        Δ = δₖ₀ - δₖ₁

        if ctr > 10 && δₖ₁ >= δₖ₀
            δₖ₀ < δₖ₁ && @warn "Increasing Δ, terminating search"
            cvg = true
            break
        end

        zₖs .-= Δₖ

        if δₖ₁^2 <= (δₖ₀ - δₖ₁) * τ
            cvg = true
            break
        end

        δₖ₀ = δₖ₁
    end
    verbose && show_stats(stats(p, zₖs, ls)...)

    if cvg
        return zₖs
    else
        @info ("""
The multiplicity count may be in error: the initial guess for the roots failed
to converge to a pejorative root.
        """)
        return(zₖs)
    end

end

# Helper: compute Π (x-zᵢ)^lᵢ directly
function evalG!(G, zs::Vector{T}, ls::Vector) where {T}

    G .= zero(T)
    G[1] = one(T)

    for (z,l) in zip(zs, ls)
        for _ in 1:l
            for j in length(G)-1:-1:1
                G[1 + j] -= z * G[j]
            end
        end
    end

    nothing
end

# For eachcolumn, computes (3.8)
# -lⱼ * (x - xⱼ)^(lⱼ-1) * Π_{k ≠ j} (x - zₖ)^lₖ
function evalJ!(J, zs::Vector{T}, ls::Vector) where {T}

    J .= 0
    n, m = sum(ls), length(zs)

    evalG!(view(J, 1:1+n-m, 1), zs, ls .- 1)

    for (j′, lⱼ) in enumerate(reverse(ls))
        for i in 1+n-m:-1:1
            J[i, m - j′ + 1] = -lⱼ * J[i, 1]
        end
    end

    for  (j, lⱼ) in enumerate(ls)
        for (k, zₖ) in enumerate(zs)
            k == j && continue
            for i in n-1:-1:1
                J[1+i,j] -= zₖ * J[i,j]
            end
        end
    end
    nothing
end

# Defn (3.5) condition number of z with respect to the multiplicity l and weight W
cond_zl(p::AbstractPolynomial, zs::Vector{S}, ls) where {S}  = cond_zl(reverse(coeffs(p)), zs, ls)

function cond_zl(p, zs::Vector{S}, ls) where {S}
    J = zeros(S, sum(ls), length(zs))
    W = diagm(0 => [min(1, 1/abs(aᵢ)) for aᵢ in p[2:end]])
    evalJ!(J, zs, ls)
    _,R = qr(W*J)
    σ = Polynomials.NGCD.smallest_singular_value(UpperTriangular(R))
    1 / σ
end

backward_error(p::AbstractPolynomial, zs::Vector{S}, ls) where {S}  =
    backward_error(reverse(coeffs(p)), zs, ls)

function backward_error(p, z̃s::Vector{S}, ls) where {S}
    # ‖ ã - a ‖w
    G = zeros(S,  1 + sum(ls))
    evalG!(G, z̃s, ls)
    a = p[2:end]./p[1]
    W = Diagonal([min(1, 1/abs(aᵢ)) for aᵢ in a])
    u = G[2:end] .- a
    norm(W*u,2)
end

function stats(p, zs, ls)
    cond_zl(p, zs, ls), backward_error(p, zs, ls)
end

function show_stats(κ, ϵ)
    println("")
    println("Condition number κ(z; l) = ", κ)
    println("Backward error ‖ã - a‖w = ", ϵ)
    println("Estimated forward root error ‖z̃ - z‖₂ = ", κ * ϵ)
    println("")
end

end
