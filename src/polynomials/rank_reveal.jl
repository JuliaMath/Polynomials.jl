# implement rank_reveal algorithm others of Zeng
# http://archive.ymsc.tsinghua.edu.cn/pacm_download/285/8743-RankRev_paper.pdf
export rank_reveal, rank_reveal!, smallest_singular_value!

using LinearAlgebra
using Random

# Numerical GCD problem
# (p,q) are a possible pertubation of (p̃, q̃) which have gcd, ũ, of degree k. Goal is to find u with min||u-kũ|| = O(||(p,q) - (p̃, q̃)||
# * numerical gcd is exact gcd of nearby problem (p̃, q̃) within ϵ
# * this nearby problem resides in the largest manifold w/in ϵ
# * for this manifold, numerical GCD is minimal distance to (p,q)

"""
    ngcd(p::PnPolynomial{T,X}, q::PnPolynomial{T,X}, [k::Int];
             ϵₐ = eps(real(T))^(5/6),
             ϵᵣ = eps(real(T)),
             ρₐ = ϵₐ,
             ρᵣ = ϵᵣ,
             λ,
             scale::Bool=false,


Return ``u, v, w, Θ, κ`` where ``u⋅v ≈ p`` and ``u⋅w ≈ q`` (polynomial
multiplication); ``Θ`` (`\\Theta[tab]`) is the residual error (``‖
[u⋅v,u⋅w] - [p,q] ‖``); and ``κ`` (`\\kappa[tab]`) is the numerical gcd
condition number estimate. When `scale=true`, ``u⋅v ≈ ps/‖ps‖₂`` and
``u⋅w ≈ qs/‖qs‖₂``.

The numerical GCD problem is defined in [1] (5.4). Let ``(p,q)`` be a
polynomial pair with degree ``m``, ``n``. Let ``Ρ_{mn}`` be set of all
such pairs. Any given pair of polynomials has an exact greatest common
divisor, ``u``, of degree ``k``, defined up to constant factors. Let
``Ρᵏmn`` be the manifold of all such ``(p,q)`` pairs with exact gcd of
degree ``k``. A given pair ``(p,q)`` with exact gcd of degree ``j``
will have some distance ``Θᵏ`` from ``Pᵏ``.  For a given threshold ``ϵ
> 0`` a numerical GCD of ``(p,q)`` within ``ϵ`` is an exact GCD of a
pair ``(p̂,q̂)`` in ``Ρᵏ`` with

``‖ (p,q) - (p̂,q̂) ‖ <= Θᵏ``, where ``k`` is the largest value for which ``Θᵏ < ϵ``.

(In the ``ϵ → 0`` limit, this would be the exact GCD.)


Suppose ``(p,q)`` is an ``ϵ`` pertubation from ``(p̂,q̂)`` where ``(p̂,q̂)`` has an exact gcd of degree ``k``, then ``degree(gcdₑ(p,q)) = k``; as ``ϵ → 0``, ``gcdₑ(p,q) → gcd(p̂, q̂)``, and

``\\limsup_{(p,q)→(p̂,q̂)} \\inf{ ‖ (u,v,w) - (û,v̂,ŵ) ‖} / ‖ (p,q) - (p̂,q̂) ‖ < κₑ(p,q)``.

``κ`` is called the numerical GCD condition number.


The Zeng algorithm proposes a degree for ``u`` and then *if* a triple
``(u,v,w)`` with ``u`` of degree ``k`` and ``(u⋅v, u⋅w)`` in ``Ρᵏmn``
can be found satisfying ``‖ (u⋅v, u⋅w) - (p,q) ‖ < ϵ`` then
``(u,v,w)`` is returned; otherwise the proposed degree is reduced and
the process repeats. If not terminated, at degree ``0`` a constant gcd
is returned.

The initial proposed degree is the first ``j``, ``j=min(m,n):-1:1``,
where ``Sⱼ`` is believed to have a singular value of ``0`` (``Sⱼ``
being related to the Sylvester matrix of `p` and `q`). The
verification of the proposed degree is done using a Gauss-Newton
iteration scheme holding the degree of ``u`` constant.

## Scaling:

If `scale=true` the gcd of `p/‖p‖₂` and `q/‖q₂` is identified. When a
polynomials has a large norm, scaling -- or using a scaled tolerance
-- can be necessary to numerically identify the degree of the gcd.

## Tolerances:

There are two places where tolerances are utilized:

* The value `ϵ` is used to determine if the Sylvester matrix, `Sⱼ` is singular. The theory states that if `Θₖ < ϵ` then `σ₋₁ < ϵ √(m - j + 1)`. We take `ϵ = max(ϵₐ, ‖p‖₂*ϵᵣ)`.

The default choice for `ϵ` works reasonably well for a range of polynomials, but scaling or some other choice of `ϵ` is needed for some cases.

* to test if ``(u ⋅ v, u ⋅ w) ≈ (p,q)`` a tolerance of `ρ` is used. The theory has `ρ = ϵ`, but this implementation allows `ρ` to be set using a value of `ρ=max(ρₐ, (‖p‖₂+‖q‖₂)*κ*ρᵣ)`, which seems to track the scaling that needed due to floating point approximations.


## Specified degree:

When `k` is specified, a value for ``(u,v,w)`` is identified with ``degree(u)=k``. No tolerances are utilized in computing ``Θᵏ``.


Output:

The function outputs a named tuple with names (`u`, `v`, `w`, `Θ`, `κ`). The components `u`,`v`,`w` estimate the gcd and give the divisors. The value is the residual error and `κ` estimates the numerical condition number.

Example:

```
using Polynomials
x = variable(Polynomial{Float64})
p = (x+10)*(x^9 + x^8/3 + 1)
q = (x+10)*(x^9 + x^8/7 - 6/7)
gcd(p, q) # u a constant
gcd(p,q, method=:numerical)  # u a degree 1 polynomial
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
              ϵₐ = eps(real(T))^(5/6),
              ϵᵣ = eps(real(T)), # backward error tolerance of (p,q), (p̃, q̃)
              ρₐ = ϵₐ,                  # residual over Πₖ
              ρᵣ = ϵᵣ,
              scale::Bool=false,
              λ = one(T),
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


    # pre-allocate storage
    Q = zeros(T, m + n, m + n) # for QR decomposition of Sylvester matrices
    R = zeros(T, m + n, m + n)
    uv = copy(p) # storage for uᵢ * vᵢ
    uw = copy(q) # storage for uᵢ * wᵢ
    x = Vector{T}(undef, m + n) # used to find σ₋₁

    # j is degree of proposed gcd j ≤ n ≤ m
    j = n  # We count down Sn, S_{n-1}, ..., S₂, S₁
    Sₓ = SylvesterMatrix(p, q, j)    # initial Sylvester matrix

    A0 = zeros(T, m+1, 2) # storage for use with extend_QR!
    A0[:,1] = coeffs(p)
    A0[end-n:end,2] = coeffs(q)

    nr, nc = size(Sₓ) # m+1, m-n+2
    F = qr(Sₓ)
    Q[1:nr, 1:nr] .= F.Q
    R[1:nc, 1:nc] .= F.R

    # tolerances
    ϵₘ = eps(real(T)) # machine tolerance
    ϵₐ, ρₐ = λ*ϵₐ, λ*ρₐ
    ϵ = max(ϵₐ, (np₂+nq₂) * ϵᵣ)


    while true

        V = UpperTriangular(view(R, 1:nc, 1:nc))
        xx = view(x, 1:nc)
        λ = norm(Q,Inf)*norm(R,Inf)
        σ₋₁ = smallest_singular_value!(xx, V, ϵ *  sqrt(m - j + 1), λ*ϵₘ)
        @show j, σ₋₁, ϵ *  sqrt(m - j + 1)

        # Lemma 7.1: If (p,q) is w/in ϵ of P^k_{mn} then σ₋₁ < ϵ√(m-j+1)
        if σ₋₁ ≤ ϵ *  sqrt(m - j + 1)
            # candidate for degree; refine u₀, vₒ, w₀ to see if ρ < ϵ
            if iszero(σ₋₁)
                @info "Determinant is zero, which shouldn't be the case. Treat results with scrutiny"
                u, v, w = initial_uvw(Val(:iszero), j, p, q, xx)
            else
                u, v, w = initial_uvw(Val(:ispossible), j, p, q, xx)
            end
            @show residual_error(p,q,u*v, u*w)
            ρₖ, κ = refine_uvw!(u, v, w, p, q, uv, uw)
            ρ = max(ρₐ, (np₂ + nq₂) * κ * ρᵣ)
            @show ρₖ, ρ, κ
            if ρₖ ≤ ρ
                return (u=u, v=v, w=w, Θ=ρₖ, κ=κ)
            end

        end

        # reduce possible degree of u and try again with Sⱼ₋₁
        # unless we hit specified minimum, in which case return it
        # minⱼ = -1
        if j == minⱼ
            u, v, w = initial_uvw(Val(:ispossible), j, p, q, xx)
            ρₖ, κ = refine_uvw!(u,v,w, p, q, uv, uw)
            return (u=u, v=v, w=w, Θ=ρₖ, κ=κ)
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
    #u, v, w = initial_uvw(Val(:ispossible), k, p, q, x)
    @show residual_error(p,q,u*v, u*w)
    ρₖ, κ = refine_uvw!(u, v, w, p, q, u*v, u*w)
    return (u=u, v=v, w=w, Θ=ρₖ, κ=κ)
end

function trivial_gcd(p::P, q) where {T, X, P <: PnPolynomial{T, X}}
    u, v, w = one(P), p, q
    return (u=u, v=v, w=w, Θ=zero(T), κ=typemax(real(T)))
end


# A = QR solve by iteration for largest eigenvalue of A^-1
# modifies w in place

# https://arxiv.org/abs/2103.04196
# Find smallest singular value
# stop when values less then atol or Δ=sⱼ - s₋₁ is small
function smallest_singular_value!(w, R::UpperTriangular{T},
                                  atol, # numerical tolerance
                                  rtol# ||A||_oo * eps(T)
                                  ) where {T}

    # Cant' handle singular matrices
    iszero(det(R)) && return zero(T)

    MAXSTEPS = 50

    sⱼ = sⱼ₋₁ = typemax(real(T))

    rand!(w)
    w ./= norm(w)

    j = 1
    while true
        sⱼ₋₁ = sⱼ
        sⱼ = smallest_singular_value_one_step!(w, R)
        δ = max(atol, norm(R, Inf)*rtol)
        if sⱼ ≤ δ || abs(sⱼ - sⱼ₋₁) ≤ sⱼ * rtol
        #if sⱼ ≤ atol || abs(sⱼ - sⱼ₋₁) ≤ sⱼ * rtol
            break
        end

        j += 1
        j >= MAXSTEPS && break
    end

    return sⱼ

end

# no tolerance; stop when improvment stops
function smallest_singular_value(A)
    R = UpperTriangular(qr(A).R)
    w = Vector{eltype(R)}(undef, size(R, 2))
    smallest_singular_value!(w, R)
end

function smallest_singular_value!(w, R::UpperTriangular{T}) where {T}
    iszero(det(R)) && return zero(T)

    MAXSTEPS = 50

    sⱼ = typemax(real(T))

    rand!(w)
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
    @show j
    m, n = length(p)-1, length(q)-1
    vᵢ = vcat(2:m-n+2, m-n+4:2:length(x))
    wᵢ = m-n+3 > length(x) ? [1] : vcat(1, (m-n+3):2:length(x))

    v = P(-x[vᵢ])
    w = P(x[wᵢ])

    # p194 3.9 C_k(v) u = p or Ck(w) u = q; this uses 10.2
    u = solve_u(v, w, p, q, j)
    return u, v, w

end

# find u₎, v₀. w₀ when R is singular
function initial_uvw(::Val{:iszero}, j, p::P, q::Q, x) where {T,X,
                                                              P<:PnPolynomial{T,X},
                                                              Q<:PnPolynomial{T,X}}

    m,n = length(p)-1, length(q)-1
    S = [convmtx(p, n-j+1) convmtx(q, m-j+1)]

    F = qr(S)
    R = UpperTriangular(F.R)

    if iszero(det(R))
        x .= eigvals(R)[:,1]
    else
        x .= ones(T, size(R,2))
        ldiv!(R', x)
        x ./= norm(x,2)
        ldiv!(R, x)
        x ./= norm(x)
    end

    w = P(x[1:n-j+1])
    v = P(-x[(n-j+2):end])

    u = solve_u(v,w,p,q,j)
    return u,v,w
end

function initial_uvw(::Val{:constant}, j, p::P, q, x) where {T,X,P<:PnPolynomial{T,X}}
    return one(P), p, q
end

# extend QR to next size
# Q gets a 1 in nc,nc, 0s should be elswhere
function extend_QR!(Q, R, nr, nc, A0)

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

    # storage
    h, β =  u, dot(u,u)  # h = constant * u₀ is used
    A = JF(h, u, v, w)
    Δfβ = Fmp(dot(h, u) - β, p, q, uv, uw)
    Δz = zeros(T, length(u) + length(v) + length(w))
    n = size(A, 2)
    R = UpperTriangular(Matrix{T}(undef, n, n))
    ũ, ṽ, w̃ = copy(u), copy(v), copy(w)

    steps = 0
    @show steps, ρ₁

    minᵢ, Maxᵢ = 3, 20



    while ρ₁ > 0.0
        steps += 1
        refine_uvw_step!(ũ, ṽ, w̃,
                         Δz, A, Δfβ, R)

        mul!(uv, ũ, ṽ)
        mul!(uw, ũ, w̃)
        ρ′ = residual_error(p, q, uv, uw)
        @show steps, ρ′
        # don't worry as much about first few,
        # but afterwards each step must be productive
        # terminate when no longer decreasing
        #if steps < minᵢ || (steps ≤ Maxᵢ && ρ′ < 0.95*ρ₁)
        if ρ′ < ρ₁ || (steps ≤ minᵢ && ρ′ ≤ 1.1*ρ₁)
            ρ₁ = ρ′
            copy!(u.coeffs, ũ.coeffs)
            copy!(v.coeffs, ṽ.coeffs)
            copy!(w.coeffs, w̃.coeffs)
            steps ≥ Maxᵢ && break
            # update A,b for next iteration
            JF!(A, h, u, v, w)
            Fmp!(Δfβ,  dot(h, u) - β, p, q, uv, uw)
        else

            break
        end
    end
    σ₂ = smallest_singular_value_one_step!(Δz, R)
    κ = 1 / σ₂
    return ρ₁, κ

end

# update u,v,w, uv, uw
# computes update step of zₖ₊₁ = zₖ - Δz; Δz = J(zₖ)⁺(fₕ(u,v,w) - [β;p;q])
function refine_uvw_step!(u, v, w,
                          Δz, J⁺, Δfβ, R)

    qrsolve!(Δz, J⁺, Δfβ, R)
    #Δz .= J⁺ \ Δfβ

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

function qrsolve!(y::Vector{T}, A, b, R) where {T}
    q = qr(A)
    R .= UpperTriangular(q.R)
    y .= q \ b
end

# Fast least-squares solver for full column rank Hessenberg-like matrices
# By Andreas Varga
function Xqrsolve!(y::Vector{T}, A, b, R) where {T <: Float64}
    # half the time of ldiv!(y, qr(A), b)
    Base.require_one_based_indexing(A)
    m, n = size(A)
    m < n && error("Column dimension exceeds row dimension")
    _, τ = LinearAlgebra.LAPACK.geqrf!(A)
    tran = T <: Complex ? 'C' : 'T'
    LinearAlgebra.LAPACK.ormqr!('L', tran, A, τ, view(b,:,1:1))
    R .= UpperTriangular(triu(A[1:n,:]))
    ldiv!(y, R, view(b,1:n))
end

function qrsolve!(y::Vector{T}, A, b) where {T <: Float64}
    n = size(A, 2)
    R = UpperTriangular(Matrix{Float64}(undef, n, n))
    qrsolve!(y, A, b, R)
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
    Sₓ = hcat(convmtx(p, n-j + 1 ),  convmtx(q, m-j + 1))
end


## ----
"""

This is the high rank-revealing algorithm
"""
function rank_reveal(A; kwargs...)
    m, n = size(A)
    #@assert m > n
    q = qr(A)

    R = UpperTriangular(q.R)
    w = Vector{eltype(A)}(undef, n)
    λ = norm(A, Inf)
    τ = norm(R, Inf)
    rank_reveal!(R, w, λ, τ; kwargs...)

end

# non allocating numeric-rank reveal
function rank_reveal!(R::LinearAlgebra.UpperTriangular{T, Matrix{T}}, w,
                      λ, τ;
                      θ = sqrt(eps(real(T))),
                      ϵ = eps(real(T))) where {T}

    n = size(R, 2)

    d = prod(R[i,i] for i ∈ 1:n)
    iszero(d) && return(0, zero(T)) # Cant' handle singular matrices


    ϵₘ = λ * ϵ
    MAXSTEPS = 50



    r = n
    sⱼ = sⱼ₋₁ = typemax(real(T))

    for k ∈ 1:n
        sⱼ = smallest_singular_value!(w, R, θ, ϵ)

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


## Tests
function test1()
    x = variable(PnPolynomial)
    α(j, n) = cospi(j/n)
    β(j, n) = sinpi(j/n)
    k(n) = n ÷ 2
    r1, r2 = 1/2, 3/2
    u(n) = prod((x - r1*α(j,n))^2 + (r1*β(j,n))^2 for j ∈ 1:k(n))
    v(n) = prod((x - r2*α(j,n))^2 + (r2*β(j,n))^2 for j ∈ 1:k(n))
    w(n) = prod((x - r1*α(j,n))^2 + (r1*β(j,n))^2 for j ∈ (k(n)+1):n)


    for n ∈ (6, 10, 16, 18, 20)
        p = u(n)*v(n)
        q = u(n)*w(n)
        a = ngcd(p, q)
        @show n, a.κ, a.Θ
    end
end


function test2()
    x = variable(PnPolynomial)
    xⱼ(j) = (-1)^j * (j/2)
    p = prod(x - xⱼ(j) for j ∈ 1:10)
    q = prod(x - xⱼ(j) + 1/10^j for j ∈ 1:10)

    for k = 2:10
        ϵ = 1/10^k
        a = ngcd(p, q, ϵₐ = ϵ)#, ϵᵣ=0.0)
        @show k, ϵ, degree(a.u), a.Θ
    end
end

function test3()
    x = variable(PnPolynomial)
    v = 1 + x + x^2 + x^3
    w = 1 - x + x^2 - x^3 + x^4


    for n ∈ (5, 20, 50, 80, 100, 200, 500, 1000, 2000)
        us =  5*(2rand(n) .- 1)
        u = PnPolynomial(us)
        p, q = u*w, u*v

        a = ngcd(p, q)
        err = norm(a.u/a.u[end] - u/u[end], Inf)
        @show n, degree(a.u), err
    end
end

function test5()
    x = variable(PnPolynomial)
    v = 1 + x + x^2 + x^3
    w = 1 - x + x^2 - x^3 + x^4

    N = 100
    errs = zeros(N)
    for i ∈ 1:N

        u = zero(x)
        for j ∈ 0:15
            cⱼ = rand(-5:5)
            iszero(cⱼ) && continue
            eⱼ = rand(0:6)
            u += cⱼ * 10^eⱼ*x^j
        end
        a = ngcd(u*w, u*v)
        û = u/u[end]
        ũ = a.u / a.u[end]

        errs[i] = norm(û - ũ, Inf)
    end

    round.(-log10.(errs), digits=1)
end

function test6()
    x = variable(PnPolynomial)
    for m ∈ ((2,1,1,0),
             (3,2,1,0),
             (4,3,2,1),
             (5,3,2,1),
             (9,6,4,2),
             (20, 14, 10, 5),
             (80,60,40,20),
             (100, 60, 40, 20))

        p = prod((x - i)^mᵢ for (i,mᵢ) ∈ enumerate(m))
        u = prod((x - i)^(max(0,mᵢ-1)) for (i,mᵢ) ∈ enumerate(m))
        a = ngcd(p, derivative(p); scale=true)

        u = u /u[end]
        ũ = a.u/a.u[end]
        Δ = norm(u - ũ, Inf)
        @show Δ, degree(a.v), a.Θ

    end
end


#==
function ngcd(p::P, q::Q,
              args...; kwargs...) where {T,X,P<:Polynomials.StandardBasisPolynomial{T,X},
                                         S,Y,Q<:Polynomials.StandardBasisPolynomial{S,Y}}

    if (degree(q) > degree(p))
        return ngcd(q, p,args...;kwargs...)
    end
    # if degree(p) > 5*(1+degree(q))
    #     a,b = divrem(p,q)
    #     return ngcd(q,b, args...; λ=100, kwargs...)
    # end

    # easy cases
    degree(p) < 0  && return (u=q,      v=p, w=one(q),  θ=NaN, κ=NaN)
    degree(p) == 0 && return (u=one(q), v=p, w=q,       θ=NaN, κ=NaN)
    degree(q) < 0  && return (u=one(q), v=p, w=zero(q), θ=NaN, κ=NaN)
    degree(q) == 0 && return (u=one(p), v=p, w=q,       Θ=NaN, κ=NaN)
    p ≈ q          && return (u=p,v=one(p),  w=one(p),  θ=NaN, κ=NaN)
    Polynomials.assert_same_variable(p,q)

    R = promote_type(float(T), float(S))
    𝑷 = Polynomials.constructorof(promote_type(P,Q)){R,X}

    ps = R[pᵢ for pᵢ ∈ coeffs(p)]
    qs = R[qᵢ for qᵢ ∈ coeffs(q)]

    # cancel zeros
    nz = min(findfirst(!iszero, ps), findfirst(!iszero, qs))
    if nz == length(qs)
        x = variable(p)
        u = x^(nz-1)
        v,w = 𝑷(ps[nz:end]), 𝑷(qs[nz:end])
        return (u=u, v=v, w=w, Θ=NaN, κ=NaN)
    end

    ## call ngcd
    p′ = PnPolynomial{R,X}(ps[nz:end])
    q′ = PnPolynomial{R,X}(qs[nz:end])
    out = ngcd(p′, q′, args...; kwargs...)

    𝑷 = Polynomials.constructorof(promote_type(P,Q)){R,X}
    u,v,w = convert.(𝑷, (out.u, out.v, out.w))
    x = variable(u)
    for _ ∈ 1:nz
        u *= x
    end
    (u=u, v=v, w=w, Θ=out.Θ, κ = out.κ)

end

# We work with PnPolynomials over AbstractFloat coefficients
==#
