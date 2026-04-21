"""
    ngcd(p, q, [k]; kwargs...)

Find numerical GCD of polynomials `p` and `q`. Refer to `NGCD.ngcd` for details.

The main entry point for this function is `gcd(p, q, method=:numerical)`, but `ngcd` outputs the gcd factorization -- `u, v, w` with `u*v ≈ p` and `u*w ≈ q` -- along with `Θ`, an estimate on how close `p`,`q` is to a gcd factorization of degree `k` and `κ` the GCD condition number.

In the case `degree(p) ≫ degree(q)`,  a heuristic is employed to first call one step of the Euclidean gcd approach, and then call `ngcd` with relaxed tolerances.

"""
function ngcd(p::P, q::Q,
              args...;
              kwargs...) where {T,X,P<:StandardBasisPolynomial{T,X},
                                S,Y,Q<:StandardBasisPolynomial{S,Y}}


    # easy cases
    degree(p) < 0  && return (u=q,      v=p, w=one(q),  θ=NaN, κ=NaN)
    degree(p) == 0 && return (u=one(q), v=p, w=q,       θ=NaN, κ=NaN)
    degree(q) < 0  && return (u=one(q), v=p, w=zero(q), θ=NaN, κ=NaN)
    degree(q) == 0 && return (u=one(p), v=p, w=q,       Θ=NaN, κ=NaN)


    if (degree(q) > degree(p))
        u,w,v,Θ,κ =  ngcd(q,p,args...; kwargs...)
        return (u=u,v=v,w=w, Θ=Θ, κ=κ)
    end

    R = promote_type(float(T))
    𝑷 = Polynomials.constructorof(P){R,X}

    ps = R[pᵢ for pᵢ ∈ coeffs(p)]
    qs = R[qᵢ for qᵢ ∈ coeffs(q)]

    # cancel zeros
    nz = min(findfirst(!iszero, ps), findfirst(!iszero, qs))
    if nz == length(qs)
        x = variable(p)
        u = x^(nz-1)
        ps,qs = ps[nz:end], qs[nz:end]
        v,w = 𝑷(ps), 𝑷(qs)
        return (u=u, v=v, w=w, Θ=NaN, κ=NaN)
    elseif 1 <= nz < length(qs)
        ps,qs = ps[nz:end], qs[nz:end]
        p,q = P(ps), Q(qs)
    end

    if degree(p) > 5*(1+degree(q))
        a,b = divrem(p,q)
        return ngcd(q, b, args...; λ=100,  kwargs...)
    end

    # other easy cases
    p ≈ q          && return (u=p,v=one(p),  w=one(p),  θ=NaN, κ=NaN)
    Polynomials.assert_same_variable(p,q)

    ## call ngcd
    P′ = PnPolynomial
    p′ = P′{R,X}(ps)
    q′ = P′{R,X}(qs)

    out = NGCD.ngcd(p′, q′, args...; kwargs...)
    ## convert to original polynomial type
    𝑷 = Polynomials.constructorof(P){R,X}
    u,v,w = convert.(𝑷, (out.u,out.v,out.w))
    if nz > 1
        u *= variable(u)^(nz-1)
    end

    (u = u, v = v, w = w, Θ = out.Θ, κ = out.κ)

end

"""
    square_free(p)

Use `ngcd` to identify the square-free part of the polynomial `p`.
"""
square_free(p) = ngcd(p, derivative(p)).v

"""
    rank_reveal(A; atol, rtol)

This is the high rank-revealing algorithm of [Lee, Li, and Zeng](http://archive.ymsc.tsinghua.edu.cn/pacm_download/285/8743-RankRev_paper.pdf) DOI: DOI 10.1007/s11075-017-0328-7.
"""
function rank_reveal(A::AbstractMatrix{T}; kwargs...) where {T <: AbstractFloat}
    m, n = size(A)
    @assert m ≥ n
    q = qr(A)

    R = UpperTriangular(q.R)
    w = Vector{eltype(A)}(undef, n)
    λ = norm(A, Inf)
    τ = norm(R, Inf)
    NGCD.rank_reveal!(R, w, λ, τ; kwargs...)

end



## ---- the work is done in this module

module NGCD
using Polynomials, LinearAlgebra
import Polynomials: constructorof, PnPolynomial


"""
    ngcd(p::PnPolynomial{T,X}, q::PnPolynomial{T,X}, [k::Int];
             atol = eps(real(T))^(5/6),       # residual over Πₖ
             rtol = eps(real(T)),
             satol = atol,
             srtol = rtol,
             λ=one(real(T)),
             scale::Bool=false
         )

Computes numerical GCD of polynomials `p` and `q`.

Returns ``u, v, w, Θ, κ`` where ``u⋅v ≈ p`` and ``u⋅w ≈ q`` (polynomial
multiplication); ``Θ`` (`\\Theta[tab]`) is the residual error (``‖
[u⋅v,u⋅w] - [p,q] ‖``); and ``κ`` (`\\kappa[tab]`) is the numerical gcd
condition number estimate. When `scale=true`, ``u⋅v ≈ ps/‖ps‖₂`` and
``u⋅w ≈ qs/‖qs‖₂``.

The numerical GCD problem is defined in [1] (5.4). Let ``(p,q)`` be a
polynomial pair with degree ``m``, ``n``. Let ``Ρₘₙ`` be set of all
such pairs. Any given pair of polynomials has an exact greatest common
divisor, ``u``, of degree ``k``, defined up to constant factors. Let
``Ρᵏₘₙ`` be the manifold of all such ``(p,q)`` pairs with exact gcd of
degree ``k``. A given pair ``(p,q)`` with exact gcd of degree ``j``
will have some distance ``Θᵏ`` from ``Pᵏ``.  For a given threshold
``ϵ>0`` a numerical GCD of ``(p,q)`` within ``ϵ`` is an exact GCD of a
pair ``(p̂,q̂)`` in ``Ρᵏ`` with

``‖ (p,q) - (p̂,q̂) ‖ ≤ Θᵏ``, where ``k`` is the largest value for
which ``Θᵏ < ϵ``.

(In the ``ϵ → 0`` limit, this would be the exact GCD.)


Suppose ``(p,q)`` is an ``ϵ`` perturbation from ``(p̂,q̂)`` where ``(p̂,q̂)`` has an exact gcd of degree ``k``, then ``degree(gcdₑ(p,q)) = k``; as ``ϵ → 0``, ``gcdₑ(p,q) → gcd(p̂, q̂)``, and

``\\limsup_{(p,q)→(p̂,q̂)} \\inf{ ‖ (u,v,w) - (û,v̂,ŵ) ‖} / ‖ (p,q) - (p̂,q̂) ‖ < κₑ(p,q)``.

``κ`` is called the numerical GCD condition number.


The Zeng algorithm proposes a degree for ``u`` and then *if* a triple
``(u,v,w)`` with ``u`` of degree ``k`` and ``(u⋅v, u⋅w)`` in ``Ρᵏₘₙ``
can be found satisfying ``‖ (u⋅v, u⋅w) - (p,q) ‖ < ϵ`` then
``(u,v,w)`` is returned; otherwise the proposed degree is reduced and
the process repeats. If not terminated, at degree ``0`` a constant gcd
is returned.

The initial proposed degree is the first ``j``, `j=min(m,n):-1:1`,
where ``Sⱼ`` is believed to have a singular value of ``0`` (``Sⱼ``
being related to the Sylvester matrix of `p` and `q`). The
verification of the proposed degree is done using a Gauss-Newton
iteration scheme holding the degree of ``u`` constant.

## Scaling:

If `scale=true` the gcd of ``p/‖p‖₂`` and ``q/‖q₂`` is identified. When the
polynomials have large norms, scaling -- or using a relative tolerance
-- can be necessary to numerically identify the degree of the gcd.

## Tolerances:

There are two places where tolerances are utilized:

* For a given `k`, the algorithm refines values `u,v,w`. The value `Θᵏ` is estimated by the difference between ``(u ⋅ v, u ⋅ w)`` and ``(p,q)``. A tolerance of `ρ` is used to test if this is smaller than specified. The arguments `atol` and `rtol` are used to compute `ϵ=max(atol, (‖(p,q)‖₂)*κ*rtol)`

* The value `ϵ` is also used to determine if the Sylvester matrix for a given `j`, `Sⱼ`, is singular. The theory has ``ϵ`` the same a above, but we this implementation uses `ρ = max(satol, ‖(p,q)‖₂*srtol)`, which seems to track the scaling that is needed due to floating point approximations. The theory states that if `Θᵏ < ϵ` then `σ₋₁ < ϵ √(m - j + 1)`.

The default choice for `ϵ` works reasonably well for a range of polynomials, but scaling or some other choice of `ϵ` is needed for some cases.

## Specified degree:

When `k` is specified, a value for ``(u, v, w)`` is identified with ``degree(u)=k``. No tolerances are utilized in computing ``Θᵏ``.


Output:

The function outputs a named tuple with names (`u`, `v`, `w`, `Θ`, `κ`). The components `u`,`v`,`w` estimate the gcd and give the divisors. The value `Θ` is the residual error and `κ` estimates the numerical condition number.

Example:

```jldoctest ngcd
julia> using Polynomials

julia> x = variable(Polynomial{Float64})
Polynomial(1.0*x)

julia> p = (x+10)*(x^9 + x^8/3 + 1);

julia> q = (x+10)*(x^9 + x^8/7 - 6/7);

julia> degree(gcd(p, q))
0

julia> degree(gcd(p, q, method=:numerical))  # u a degree 1 polynomial
1
```

This example perturbs `q` more than the default tolerance, so `atol` is set. We can see that the residual error found is on the order of `1e-9`:

```jldoctest ngcd
julia> p = (x-1)*(x-2)*(x-3);

julia> q = (x-2)*(x-3)*(x-4) + (1e-8)*x^2;

julia> out = Polynomials.ngcd(p, q, atol=1e-8);

julia> degree(out.u)
2

julia> round(log10(out.Θ))
-9.0

julia> out = Polynomials.ngcd(p, q);

julia> degree(out.u)
0
```


Reference:

[1] The Numerical Greatest Common Divisor of Univariate Polynomials
by Zhonggang Zeng;
[url](http://homepages.neiu.edu/~zzeng/uvgcd.pdf);
[doi](https://doi.org/10.1090/conm/556/11014)

Note: Based on work by Andreas Varga

"""
function ngcd(p::PnPolynomial{T,X},
              q::PnPolynomial{T,X};
              atol = eps(real(T))^(5/6),       # residual over Θᵏ
              rtol = Base.rtoldefault(real(T)),
              satol = atol,                    # singular tolerance
              srtol =  eps(real(T)),
              scale::Bool=false,
              λ::Real = one(real(T)),          # not used
              minⱼ = -1
              ) where {T, X}

    m,n = length(p)-1, length(q)-1
    (m == 1 || n == 0) && return trivial_gcd(p, q)

    @assert m >= n

    # scale
    np₂, nq₂ = norm(p,2), norm(q,2)
    if scale
        p ./= np₂
        q ./= nq₂
        np₂ =  nq₂ = one(T)
    end
    npq₂ = sqrt(np₂^2 + nq₂^2)


    # pre-allocate storage
    Q = zeros(T, m + n, m + n) # for QR decomposition of Sylvester matrices
    R = zeros(T, m + n, m + n)
    uv = copy(p) # storage for uᵢ * vᵢ
    uw = copy(q) # storage for uᵢ * wᵢ
    x = Vector{T}(undef, m + n) # used to find σ₋₁

    # j is degree of proposed gcd j ≤ n ≤ m
    j = n  # We count down Sn, S_{n-1}, ..., S₂, S₁
    Sₓ = SylvesterMatrix(p, q, j)    # initial Sylvester matrix [Cₙ₋ⱼ₊₁(p), Cₘ₋ⱼ₊₁(q)]

    A0 = zeros(T, m+1, 2) # storage for use with extend_QR!
    A0[:,1] = coeffs(p)
    A0[end-n:end,2] = coeffs(q)

    nr, nc = size(Sₓ) # m+1, m-n+2
    F = qr(Sₓ)
    dest = view(Q, 1:nr, 1:nr)
    copyto!(dest, I)
    lmul!(F.Q, dest)
    R[1:nc, 1:nc] .= F.R

    # tolerances
    atol, satol, rtol, srtol = λ*atol, λ*satol, λ*rtol, λ*srtol
    ρ = max(satol, npq₂ * srtol)

    while true
        V = UpperTriangular(view(R, 1:nc, 1:nc))
        xx = view(x, 1:nc)
        σ₋₁ = smallest_singular_value!(xx, V, ρ *  sqrt(m - j + 1))
        #@show j, σ₋₁, ρ *  sqrt(m - j + 1), npq₂

        # Lemma 7.1: If (p,q) is w/in ϵ of P^k_{mn} then σ₋₁ < ϵ√(m-j+1)
        if σ₋₁ ≤ ρ *  sqrt(m - j + 1)
            # candidate for degree; refine u₀, vₒ, w₀ to see if ρ < ϵ
            if iszero(σ₋₁)
                # determinant is 0
                u, v, w = initial_uvw(Val(:iszero), j, p, q, xx)
            else
                u, v, w = initial_uvw(Val(:ispossible), j, p, q, xx)
            end
            ϵₖ, κ = refine_uvw!(u, v, w, p, q, uv, uw)

            # we have limsup Θᵏ / ‖(p,q) - (p̃,q̃)‖ = κ, so
            # ‖Θᵏ‖ ≤ κ ⋅ ‖(p,q)‖ ⋅ ϵ seems a reasonable heuristic.
            # Too tight a tolerance and the right degree will be missed; too
            # lax, and larger degrees will be accepted. We are using
            # `√eps()` for `rtol`, but that may be too lax and is subject to
            # change.
            ϵ = max(atol, npq₂ * κ * rtol)
            #@show ϵₖ, ϵ, κ
            if ϵₖ ≤ ϵ
                return (u=u, v=v, w=w, Θ=ϵₖ, κ=κ)
            end
            #@show :failure, j
        end

        # reduce possible degree of u and try again with Sⱼ₋₁
        # unless we hit specified minimum, in which case return it
        # minⱼ = -1
        if j == minⱼ
            u, v, w = initial_uvw(Val(:ispossible), j, p, q, xx)
            ϵₖ, κ = refine_uvw!(u, v ,w, p, q, uv, uw)
            return (u=u, v=v, w=w, Θ=ϵₖ, κ=κ)
        end

        # Try again with a smaller j
        j -= 1
        nr += 1
        nc += 2
        nc > nr && break
        extend_QR!(Q, R, nr, nc, A0) # before Q⋅R = Sⱼ, now Q⋅R = Sⱼ₋₁
    end

    return trivial_gcd(p, q)

end

# fixed `k`
function ngcd(p::P, q::P, k::Int) where {T, X, P <: PnPolynomial{T,X}}

    m, n = length(p)-1, length(q)-1
    Sₓ = SylvesterMatrix(p,q,k)
    F = qr(Sₓ)
    R = UpperTriangular(F.R)
    x = zeros(T, size(Sₓ, 2))
    np = norm(p)
    σ₋₁ = smallest_singular_value!(x, R)
    w, v = P(x[1:(n-k+1)]), P(-x[(n-k+2):end])
    u = solve_u(v, w, p, q, k)
    ρₖ, κ = refine_uvw!(u, v, w, p, q, u*v, u*w)
    return (u=u, v=v, w=w, Θ=ρₖ, κ=κ)
end

function trivial_gcd(p::P, q) where {T, X, P <: PnPolynomial{T, X}}
    u, v, w = one(P), p, q
    return (u=u, v=v, w=w, Θ=zero(T), κ=NaN)
end


# A = QR solve by iteration for largest eigenvalue of A^-1
# modifies w in place

# https://arxiv.org/abs/2103.04196
# Find smallest singular value
# stop when value is less then θ or Δ=sⱼ - s₋₁ is small
function smallest_singular_value!(w, R::UpperTriangular{T},
                                  θ,
                                  ϵₘ = eps(real(T))
                                  ) where {T}

    # Cant' handle singular matrices
    iszero(det(R)) && return zero(T)

    nRₒₒ = norm(R, Inf)
    MAXSTEPS = 50
    sⱼ = sⱼ₋₁ = typemax(real(T))

    w .= one(T)
    w ./= norm(w)

    j = 1
    while true
        sⱼ₋₁ = sⱼ
        sⱼ = smallest_singular_value_one_step!(w, R)

        if sⱼ ≤ θ || abs(sⱼ - sⱼ₋₁) ≤ max(1, sⱼ)  * ϵₘ *  nRₒₒ
            break
        end

        j += 1
        j >= MAXSTEPS && break
    end

    return sⱼ

end

# no tolerance; stop when improvement stops
function smallest_singular_value(A)
    R = UpperTriangular(qr(A).R)
    w = Vector{eltype(R)}(undef, size(R, 2))
    smallest_singular_value!(w, R)
end

function smallest_singular_value!(w, R::UpperTriangular{T}) where {T}
    iszero(det(R)) && return zero(T)

    MAXSTEPS = 50

    sⱼ = typemax(real(T))

    w .= one(T)
    w ./= norm(w)

    j = 1

    w̃ = copy(w)
    while true
        s′ = smallest_singular_value_one_step!(w̃, R)

        s′ >  sⱼ && break

        copy!(w, w̃)
        sⱼ = s′
        j += 1
        j > MAXSTEPS && break

    end

    return sⱼ

end


# modify w, return sⱼ after one step
# uses R from QR factorization
function smallest_singular_value_one_step!(w, R)
    ldiv!(R', w)
    w ./= norm(w,2)
    ldiv!(R, w)
    sⱼ = 1/norm(w, 2)
    w .*= sⱼ
    return sⱼ
end


# solve for u from [v,w] \ [p,q]
function solve_u(v::P, w, p, q, j) where {T, X, P<:PnPolynomial{T,X}}
    A = [convmtx(v, j+1); convmtx(w, j+1)]
    b = vcat(coeffs(p), coeffs(q))
    u = A \ b
    return P(u)
end

## Find u₀,v₀,w₀ from right singular vector
function initial_uvw(::Val{:ispossible}, j, p::P, q::Q, x) where {T,X,
                                                              P<:PnPolynomial{T,X},
                                                                  Q<:PnPolynomial{T,X}}
    # Sk*[w;-v] = 0, so pick out v,w after applying permutation
    m, n = length(p)-1, length(q)-1

    vᵢ = vcat(2:m-n+2, m-n+4:2:length(x))
    wᵢ = m-n+3 > length(x) ? [1] : vcat(1, (m-n+3):2:length(x))

    v = P(-x[vᵢ])
    w = P(x[wᵢ])

    # p194 3.9 C_k(v) u = p or Ck(w) u = q; this uses 10.2
    u = solve_u(v, w, p, q, j)
    return u, v, w

end

# find u₀, v₀. w₀ when R is singular.
function initial_uvw(::Val{:iszero}, j, p::P, q::Q, x) where {T,X,
                                                              P<:PnPolynomial{T,X},
                                                              Q<:PnPolynomial{T,X}}

    m,n = length(p)-1, length(q)-1
    S = SylvesterMatrix(p,q,j)
    F = qr(S)
    R = UpperTriangular(F.R)

    if iszero(det(R))
        x .= eigvals(R)[:,1]
    else
        x .= one(T)
        smallest_singular_value!(x, R)
    end

    w = P(x[1:n-j+1]) # ordering of S is not interlaced
    v = P(-x[(n-j+2):end])

    u = solve_u(v,w,p,q,j)
    return u,v,w
end


# extend QR to next size
# Q gets a 1 in nc,nc, 0s should be elsewhere
function extend_QR!(Q::AbstractMatrix{T}, R::AbstractMatrix{T}, nr, nc, A0) where {T}

    #old Q is m x m, old R is n x n; we add to these
    n = nc - 2
    m = nr - 1
    k,l = size(A0)

    # add two columns to R
    # need to apply Q to top part of new columns
    R[nr-k+1:nr, (nc-1):nc] = A0
    R[1:nr-1, (nc-1):nc] = (view(Q, 1:nr-1, 1:nr-1))' *  R[1:nr-1, (nc-1):nc]

    # extend Q with row and column with identity
    Q[nr,nr] = 1

    # Make R upper triangular using Givens rotations
    for j in nr-1:-1:nc-1
        gj,_ = givens(R[j,nc-1], R[j+1,nc-1], j, j+1)
        rmul!(Q, gj')
        lmul!(gj, R)
    end

    for j in nr-1:-1:nc
        gj,_ = givens(R[j,nc], R[j+1,nc], j, j+1)
        rmul!(Q, gj')
        lmul!(gj, R)
    end

    return nothing

end

## attempt to refine u,v,w
## return residual error, ρ, estimate for 1/σ_2, κ
function refine_uvw!(u::P, v::P, w::P,
                     p, q, uv, uw) where {T,X,
                                          P<:PnPolynomial{T,X}}
    mul!(uv, u, v)
    mul!(uw, u, w)
    ρ₁ = residual_error(p, q, uv, uw)
    iszero(ρ₁) && return (ρ₁, NaN)
    # storage
    h, β =  u, dot(u.coeffs,u.coeffs)  # h = constant * u₀ is used
    A = JF(h, u, v, w)
    Δfβ = Fmp(dot(h.coeffs, u.coeffs) - β, p, q, uv, uw)
    Δz = ones(T, length(u) + length(v) + length(w))
    n = size(A, 2)
    R = UpperTriangular(Matrix{T}(undef, n, n))
    R′ = deepcopy(R)
    ũ, ṽ, w̃ = copy(u), copy(v), copy(w)

    steps = 0
    #@show steps, ρ₁
    minᵢ, Maxᵢ = 3, 20

    while ρ₁ > 0.0
        steps += 1
        refine_uvw_step!(ũ, ṽ, w̃,
                         Δz, A, Δfβ, R)
        mul!(uv, ũ, ṽ)
        mul!(uw, ũ, w̃)
        ρ′ = residual_error(p, q, uv, uw)
        #@show steps, ρ′
        # don't worry as much about first few,
        # but afterwards each step must be productive
        # terminate when no longer decreasing
        #if steps < minᵢ || (steps ≤ Maxᵢ && ρ′ < 0.95*ρ₁)
        if ρ′ < ρ₁ || (steps ≤ minᵢ && ρ′ ≤ 1.1*ρ₁)
            ρ₁ = ρ′
            copy!(R′, R)
            copy!(u.coeffs, ũ.coeffs)
            copy!(v.coeffs, ṽ.coeffs)
            copy!(w.coeffs, w̃.coeffs)
            steps ≥ Maxᵢ && break
            # update A,b for next iteration
            JF!(A, h, u, v, w)
            Fmp!(Δfβ,  dot(h.coeffs, u.coeffs) - β, p, q, uv, uw)
        else
            steps == 1 && copy!(R′, R)
            break
        end
    end
    smallest_singular_value_one_step!(Δz, R′) # two steps, not one
    σ₂ = smallest_singular_value_one_step!(Δz, R′)
    κ = 1/σ₂
    return ρ₁, κ

end

# update u,v,w, uv, uw
# computes update step of zₖ₊₁ = zₖ - Δz; Δz = J(zₖ)⁺(fₕ(u,v,w) - [β;p;q])
function refine_uvw_step!(u, v, w,
                          Δz, J⁺, Δfβ, R)

    qrsolve!(Δz, J⁺, Δfβ, R) # Δz .= J⁺ \ Δfβ

    m,n,l = length(u)-1, length(v)-1, length(w)-1

    Δuᵢ = 1:(m+1)
    Δvᵢ = (m+1+1):(m+1 + n+1)
    Δwᵢ = (m + 1 + n + 1 + 1):(m + n + l + 3)

    u .-= view(Δz, Δuᵢ)
    v .-= view(Δz, Δvᵢ)
    w .-= view(Δz, Δwᵢ)

end

## Jacobian of F(u,v,w) = [p, p'] is J(u,v,w)
## [h      0      0;
##  Cₖ(v) Cₘ₋ₖ(u) 0;
##  Cₖ(w)  0    Cₙ₋ₖ(u)]
function JF(h, u::P, v, w) where {T,X,P<:AbstractPolynomial{T,X}}

    # compute size needed to store
    m, k, j = length(u)-1, length(v)-1, length(w)-1
    n, l = m + k, m + j

    ai, aj = convmtx_size(v, m + 1)
    bi, bj = convmtx_size(u, k + 1)
    di, dj = convmtx_size(w, m + 1)
    fi, fj = convmtx_size(u, j + 1)
    ci, cj = ai, fj
    ei, ej = di, bj

    m, n = 1 + ai + di, aj + bj + cj

    A = zeros(T, m, n)
    JF!(A, h, u, v, w)
    A
end

# compute Jacobian of fₕ
function JF!(M, h,  u::P, v, w) where {T,X,P<:AbstractPolynomial{T,X}}

    k = length(u) - 1
    dᵥ, dᵥᵥ = length(v) - 1, length(w) - 1
    m = k + dᵥ  # degree(p)
    n = k + dᵥᵥ # degree(q)

    M[1, 1:length(h)] .= h'

    r11, c11 = convmtx_size(v, k+1) # Ck(v) size
    J11 = view(M, 2:(1+r11), 1:c11)
    convmtx!(J11, v, k+1)

    r12, c12 = convmtx_size(u, m-k+1)
    J12 = view(M, 2:(1+r11), (c11+1):(c11+1+c12))
    convmtx!(J12, u, m-k+1)

    r21, c21 = convmtx_size(w, k+1)
    J21 = view(M, (1 + r11 + 1):(1 + r11 + r21), 1:c21)
    convmtx!(J21, w, k+1)

    r23, c23 = convmtx_size(u, n - k+1)
    J23 = view(M, (1 + r11 + 1):(1 + r11 + r23),
               (c21 + c12 + 1):(c21 + c12 + c23))
    convmtx!(J23, u, n-k+1)

    return nothing
end

# create storage for fₕ(u,v,w) - [β; p; q]; fill in
function Fmp(Δ, p::PnPolynomial{T,X}, q, p̃, q̃) where {T,X}
    b = zeros(T, 1 + length(p) + length(q))
    Fmp!(b, Δ, p, q, p̃, q̃)
    b
end

## fₕ - [β; p; q]
## Δ = h⋅uᵢ - β
function Fmp!(b, Δ, p, q, uv, uw)
    b[1] = Δ
    for (i, pᵢ) ∈ pairs(p)
        b[2 + i] = uv[i] - pᵢ
    end
    for (i, qᵢ) ∈ pairs(q)
        b[2 + length(p) + i] = uw[i] - qᵢ
    end
    return nothing
end

# find ||(p,q) - (p̃, q̃)|| treating (,) as vector concatenation of coefficients
function residual_error(p::P, q, p̃, q̃) where {T,X,P<:AbstractPolynomial{T,X}}
    tot = zero(real(T))
    for (pᵢ, p̃ᵢ) in zip(p, p̃)
        tot += abs2(pᵢ - p̃ᵢ)
    end
    for (qᵢ, q̃ᵢ) in zip(q, q̃)
        tot += abs2(qᵢ - q̃ᵢ)
    end
    sqrt(tot)
end

## ---- utils
## ---- QR factorization
function qrsolve!(y::Vector{T}, A, b) where {T <: Float64}
    n = size(A, 2)
    R = UpperTriangular(Matrix{Float64}(undef, n, n))
    qrsolve!(y, A, b, R)
end

function qrsolve!(y::Vector{T}, A, b, R) where {T}
    q = qr(A)
    R .= UpperTriangular(q.R)
    y .= q \ b
end

# Fast least-squares solver for full column rank Hessenberg-like matrices
# By Andreas Varga
function qrsolve!(y::Vector{T}, A, b, R) where {T <: Float64}
    # half the time of ldiv!(y, qr(A), b)
    Base.require_one_based_indexing(A)
    m, n = size(A)
    m < n && error("Column dimension exceeds row dimension")
    _, τ = LinearAlgebra.LAPACK.geqrf!(A)
    R .= UpperTriangular(triu(A[1:n,1:n]))

    tran = T <: Complex ? 'C' : 'T'
    LinearAlgebra.LAPACK.ormqr!('L', tran, A, τ, view(b,:,1:1))
    R = UpperTriangular(triu(A[1:n,:]))
    ldiv!(y, R, view(b,1:n))
end


"""
    convmtx(v, n::Int)
    convmtx!(M,v,n)
    convmtx_size(v,n)

Convolution matrix.
C = convmtx(v,n) returns the convolution matrix C for a vector v.
If q is a column vector of length n, then C*q is the same as conv(v,q).

"""
function convmtx!(C, v::AbstractVector{T}, n::Int) where T

    #   Form C as the Toeplitz matrix
    #   C = Toeplitz([v; zeros(n-1)],[v[1]; zeros(n-1));  put Toeplitz code inline

    nv = length(v)-1

    @inbounds for j = 1:n
        C[j:j+nv,j] = v
    end

    nothing

end
convmtx_size(v::AbstractVector, n) = (n + length(v) - 1, n)
function convmtx(v::AbstractVector{T}, n::Int) where {T}
    C = zeros(T, convmtx_size(v, n)...)
    convmtx!(C, v, n)
    C
end

# multroot uses vector/matrix interface.
convmtx!(C, v::AbstractPolynomial, n::Int) = convmtx!(C, coeffs(v), n)
convmtx_size(v::AbstractPolynomial, n) = (n + length(v)-1, n)
function convmtx(v::AbstractPolynomial{T}, n::Int) where {T}
    d = length(v)-1
    C = zeros(T, (n + d, n))
    convmtx!(C, v, n)
    C
end

function SylvesterMatrix(p, q, j)
    m, n = length(p)-1, length(q) - 1
    Sₓ = hcat(convmtx(p, n - j + 1 ),  convmtx(q, m - j + 1))
end

## ----

# non allocating numeric-rank reveal
function rank_reveal!(R::LinearAlgebra.UpperTriangular{T}, w,
                      nAₒₒ, nRₒₒ;
                      atol = Base.rtoldefault(real(T)),
                      rtol = eps(real(T))) where {T <: AbstractFloat}

    n = size(R, 2)

    d = prod(R[i,i] for i ∈ 1:n)
    iszero(d) && return(0, zero(T)) # Cant' handle singular matrices


    MAXSTEPS = 50

    θ = max(atol, nAₒₒ * rtol)

    ϵₘ = nAₒₒ * eps(real(T))
    r = n
    sⱼ = sⱼ₋₁ = typemax(real(T))

    for k ∈ 1:n
        sⱼ = smallest_singular_value!(w, R, θ, ϵₘ)

        if sⱼ > θ
            break
        end

        r -= 1

        # W = [w W]
        # can eliminate allocations here!
        # RR = [τ*w'; R]
        # for i ∈ 1:n
        #     G,_ = givens(RR, i, i+1, i)
        #     lmul!(G, RR)
        # end
        # R = UpperTriangular(RR[1:end-1, :])

        τ = nRₒₒ
        a,b = τ*w[1], R[1,1]
        g,d = givens(a, b, 1, 2)
        R[1,1] = d
        for j in 2:n
            a, b =  τ*w[j], R[1, j]
            R[1,j] = conj(g.c)*a + conj(g.s)*b
            w[j] = R[2, j]
            R[2,j] = -conj(g.s)*a + conj(g.c)*b
        end

        for i ∈ 2:(n-1)
            a, b = R[i, i], w[i]
            g, d = givens(a, b, 1, 2)
            R[i,i] = d
            for j ∈ (i+1):n
                a, b =  R[i,j], w[j]
                R[i,j] = conj(g.c)*a + conj(g.s)*b
                w[j] = R[i+1, j]
                R[i+1, j] = -conj(g.s)*a + conj(g.c)*b
            end
        end

        # n
        g, d = givens(R[n,n], w[n], n-1,n)
        R[n,n] = d

    end
    r, sⱼ
end


end
