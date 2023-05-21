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

* `ρ`: the initial residual tolerance, set to `1e-13`. This is passed
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


# The multiplicity structure, `l`, gives rise to a pejorative manifold
# `Πₗ = {Gₗ(z) | z ∈ Cᵐ}`. This first finds the approximate roots, then
# finds the multiplicity structure

# compute initial root approximations with multiplicities
# without any information

# With method = :direct, use the direct method of Brehard, Poteaux, and Soudant
# based on the cofactors v,w s.t. p = u*v and q = u*w

# With method = :iterative, use the iterative strategy of Zeng
# to recover the multiplicities associated to the computed roots

# if roundmul = true, multiplicities are rounded
# to the closest integer

# Better performing :direct method by Florent Bréhard, Adrien Poteaux, and Léo Soudant [Validated root enclosures for interval polynomials with multiplicities](preprint)

function pejorative_manifold(
    p::Polynomials.StandardBasisPolynomial{T,X};
    method = :iterative,
    roundmul = true,
    θ = 1e-8,  # zero singular-value threshold
    ρ = 1e-13, #1e-10, # initial residual tolerance
    ϕ = 100,   # residual growth factor
#    cofactors=false, # return cofactors v,w?
    )  where {T,X}

    S = float(T)
    u = Polynomials.PnPolynomial{S,X}(S.(coeffs(p)))

    nu₂ = norm(u, 2)
    θ2, ρ2 =  θ * nu₂, ρ * nu₂

    u, v, w, ρⱼ, κ = Polynomials.ngcd(u, derivative(u),
        satol = θ2, srtol = zero(real(T)),
        atol = ρ2, rtol = zero(real(T)))
    ρⱼ /= nu₂

    # root approximations
    zs = roots(v)

    # recover multiplicities
    if method == :iterative
        ls = pejorative_manifold_iterative_multiplicities(
            u, zs, ρⱼ, θ, ρ, ϕ)
    elseif method == :direct
        ls = pejorative_manifold_direct_multiplicities(
            v, w, zs, roundmul=roundmul)
    else
        error("pejorative_manifold: unknown multiplicity recovering method: " * string(method))
    end

#    if cofactors
#        return zs, ls, v, w
#    else
    return zs, ls
#    end
end


# compute initial root approximations with multiplicities
# when the number k of distinct roots is known
# If method=direct and leastsquares=true, compute the cofactors v,w
# using least-squares rather than Zeng's AGCD refinement strategy
function pejorative_manifold(p::Polynomials.StandardBasisPolynomial{T,X},
    k::Int;
    method = :iterative,
    leastsquares = false,
    roundmul = true,
    θ = 1e-8,  # zero singular-value threshold
    ρ = 1e-10, # initial residual tolerance
    ϕ = 100,   # residual growth factor
#    cofactors=false, # return cofactors v,w?
    )  where {T,X}

    S = float(T)
    u = Polynomials.PnPolynomial{S,X}(S.(coeffs(p)))

    nu₂ = norm(u, 2)

    if method == :direct && leastsquares
        n = degree(u)
        Sy = Polynomials.NGCD.SylvesterMatrix(u, derivative(u), n-k)
        b = Sy[1:end-1,2*k+1] - n * Sy[1:end-1,k] # X^k*p' - n*X^{k-1}*p
        A = Sy[1:end-1,1:end .∉ [[k,2*k+1]]]
        x = zeros(S, 2*k-1)
        Polynomials.NGCD.qrsolve!(x, A, b)
        # w = n*X^{k-1} + ...
        w = Polynomials.PnPolynomial([x[1:k-1]; n])
        # v = X^k + ...
        v = Polynomials.PnPolynomial([-x[k:2*k-1]; 1])

    else
        u, v, w, ρⱼ, κ = Polynomials.ngcd(u, derivative(u),
            degree(u)-k)
        ρⱼ /= nu₂
    end

    # root approximations
    zs = roots(v)

    # recover multiplicities
    if method == :iterative
        ls = pejorative_manifold_iterative_multiplicities(
            u, zs, ρⱼ, θ, ρ, ϕ)
    elseif method == :direct
        ls = pejorative_manifold_direct_multiplicities(
            v, w, zs, roundmul=roundmul)
    else
        error("pejorative_manifold: unknown multiplicity recovering method: " * string(method))
    end

#    if cofactors
#        return zs, ls, v, w
#    else
    return zs, ls
#    end
end


# compute initial root approximations with multiplicities
# when the multiplicity structure l is known
# If method=direct and leastsquares=true, compute the cofactors v,w
# using least-squares rather than Zeng's AGCD refinement strategy
function pejorative_manifold(p::Polynomials.StandardBasisPolynomial{T,X},
    l::Vector{Int};
    method = :iterative,
    leastsquares = false,
    roundmul = true,
#    cofactors=false, # return cofactors v,w?
    )  where {T,X}

    S = float(T)
    u = Polynomials.PnPolynomial{S,X}(S.(coeffs(p)))

    # number of distinct roots
    k = sum(l .> 0)

    if method == :direct && leastsquares
        n = degree(u)
        Sy = Polynomials.NGCD.SylvesterMatrix(u, derivative(u), n-k)
        b = Sy[1:end-1,2*k+1] - n * Sy[1:end-1,k] # X^k*p' - n*X^{k-1}*p
        A = Sy[1:end-1,1:end .∉ [[k,2*k+1]]]
        x = zeros(S, 2*k-1)
        Polynomials.NGCD.qrsolve!(x, A, b)
        # w = n*X^{k-1} + ...
        w = Polynomials.PnPolynomial([x[1:k-1]; n])
        # v = X^k + ...
        v = Polynomials.PnPolynomial([-x[k:2*k-1]; 1])

    else
        u, v, w, _, _ = Polynomials.ngcd(u, derivative(u),
            degree(u)-k)
    end

    # root approximations
    zs = roots(v)

    # recover multiplicities
    if method == :iterative
        ls = pejorative_manifold_iterative_multiplicities(
            u, zs, l, roundmul=roundmul)
    elseif method == :direct
        ls = pejorative_manifold_direct_multiplicities(
            v, w, zs, l, roundmul=roundmul)
    else
        error("pejorative_manifold: unknown multiplicity recovering method: " * string(method))
    end

#    if cofactors
#        return zs, ls, v, w
#    else
    return zs, ls
#    end
end


# recover the multiplicity of each root approximation
# using the iterative method of Zeng
# but knowing the total multiplicity structure,
# hence the degree of the approximate GCD at each iteration is known
# if roundmul=true, round the floating-point multiplicities
# to the closest integer
function pejorative_manifold_iterative_multiplicities(
    u::Polynomials.PnPolynomial{T},
    zs, l::Vector{Int};
    roundmul = true) where {T}

    ls = ones(Int, length(zs))
    ll = copy(l)

    while !Polynomials.isconstant(u)
        # remove 1 to all multiplicities
        ll .-= 1
        deleteat!(ll, ll .== 0)
        k = length(ll)
        u, v, w, ρⱼ, κ = Polynomials.ngcd(u, derivative(u), degree(u)-k)
        for z ∈ roots(v)
            _, ind = findmin(abs.(zs .- z))
            ls[ind] += 1
        end
    end

    # if roundmul
    #     for i ∈ 1:length(ls)
    #         _, ind = findmin(abs.(l .- ls[i]))
    #         ls[i] = l[ind]
    #     end
    # end

    if roundmul
        ls = Int.(round.(real.(ls)))
    end

    return ls
end


# recover the multiplicity of each root approximation
# directly from the cofactors v,w s.t. p = u*v and q = u*w
# if roundmul=true, round the floating-point multiplicities
# to the closest integer
function pejorative_manifold_direct_multiplicities(
    v::Polynomials.PnPolynomial{T}, w::Polynomials.PnPolynomial{T},
    zs;
    roundmul = true) where {T}

    dv = derivative(v)
    ls = w.(zs) ./ dv.(zs)

    if roundmul
        ls = Int.(round.(real.(ls)))
    end

    return ls
end

# recover the multiplicity of each root approximation
# directly from the cofactors v,w s.t. p = u*v and q = u*w
# if roundmul=true, round the floating-point multiplicities
# to the closest integer
function pejorative_manifold_direct_multiplicities(
    v::Polynomials.PnPolynomial{T}, w::Polynomials.PnPolynomial{T},
    zs, l::Vector{Int};
    roundmul = true) where {T}

    ls = pejorative_manifold_direct_multiplicities(v, w, zs,
        roundmul = false)

    # if roundmul
    #     for i ∈ 1:length(ls)
    #         _, ind = findmin(abs.(l .- ls[i]))
    #         ls[i] = l[ind]
    #     end
    #     ls = Int.(ls)
    # end

    if roundmul
        ls = Int.(round.(real.(ls)))
    end

    return ls
end


# recover the multiplicity of each root approximation
# using the iterative method of Zeng
function pejorative_manifold_iterative_multiplicities(u::Polynomials.PnPolynomial{T},
                                            zs, ρⱼ,
                                            θ, ρ, ϕ) where {T}
    nrts = length(zs)
    ls = ones(Int, nrts)

    while !Polynomials.isconstant(u)

        nu₂ = norm(u,2)
        θ = θ * nu₂
        ρ =  max(ρ, ϕ * ρⱼ)
        ρ′ = ρ * nu₂
        u, v, w, ρⱼ, κ = Polynomials.ngcd(u, derivative(u),
                                           satol = θ, srtol = zero(real(T)),
                                           atol = ρ′, rtol = zero(real(T)),
                                           minⱼ = degree(u) - nrts
                                          )
        ρⱼ /= nu₂
        nrts = degree(v)
        for z ∈ roots(v)
            _, ind = findmin(abs.(zs .- z))
            ls[ind] += 1
        end

    end

    ls

end


"""
    pejorative_root(p, zs, ls; kwargs...)

Find a *pejorative* *root* for `p` given multiplicity structure `ls` and initial guess `zs`.

The pejorative manifold for a multiplicity structure `l` is denoted `{Gₗ(z) | z ∈ Cᵐ}`. A pejorative
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
        @info """
The multiplicity count may be in error: the initial guess for the roots failed
to converge to a pejorative root.
"""
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
    evalJ!(J, zs, ls)
    W = diagm(0 => [min(1, 1/abs(aᵢ)) for aᵢ in p[2:end]])

    σ = Polynomials.NGCD.smallest_singular_value(W*J)
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
