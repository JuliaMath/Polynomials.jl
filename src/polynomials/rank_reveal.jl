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


function Ngcd(p::PnPolynomial{T,X},
              q::PnPolynomial{T,X};
              ϵ = sqrt(eps(real(T))), # closeness of (p,q), (p̃, q̃)
              ρ = ϵ,                  # residual over Πₖ
              scale::Bool=false,
              minⱼ = -1
              ) where {T, X}


    rtol = eps(real(T))

    m,n = length(p)-1, length(q)-1
    (m == 1 || n == 0) && return trivial_gcd(p, q)

    @assert m >= n

    # scale
    if scale
        p ./= norm(p, 2)
        q ./= norm(q, 2)
    end


    # pre-allocate storage
    Q = zeros(T, m + n, m + n) # for QR decomposition of Sylvester matrices
    R = zeros(T, m + n, m + n)
    uv = copy(p) # storage for uᵢ * vᵢ
    uw = copy(q) # storage for uᵢ * wᵢ
    x = Vector{T}(undef, m + n) # used to find σ₋₁
    A0 = zeros(T, m+1, 2) # storage for use with extend_QR!
    A0[:,1] = coeffs(p)
    A0[end-n:end,2] = coeffs(q)

    # initial Sylvester matrix
    Sₓ = hcat(convmtx(p,1),  convmtx(q, m-n+1))
    nr, nc = size(Sₓ) # m+1, m-n+2
    λ = norm(Sₓ, Inf)
    F = qr(Sₓ)
    Q[1:nr, 1:nr] .= F.Q
    R[1:nc, 1:nc] .= F.R

    j = n  # We count down Sn, S_{n-1}, ..., S₂, S₁

    while true

        V = UpperTriangular(view(R, 1:nc, 1:nc))
        xx = view(x, 1:nc)

        σ₋₁ = smallest_singular_value!(xx, V, ϵ *  sqrt(m - j + 1), λ * rtol)
        @show j, σ₋₁, ϵ *  sqrt(m - j + 1), λ * rtol

        if σ₋₁ ≤ ϵ *  sqrt(m - j + 1)
            # candidate for degree; refine u₀, vₒ, w₀ to see if ρ < ϵ
            if iszero(σ₋₁)
                u, v, w = initial_uvw(Val(:iszero), j, p, q, xx)
            else
                u, v, w = initial_uvw(Val(:ispossible), j, p, q, xx)
            end
            ρₖ, κ = refine_uvw!(u, v, w, p, q, uv, uw)
            @show ρₖ ≤ ρ, ρₖ
            if ρₖ ≤ ρ
                @show :keeper, j, ρₖ
                return (u=u, v=v, w=w, θ=ρₖ, κ=κ)
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
        λ = norm(Q, Inf) * norm(R, Inf) # cheaper than ||Q*R||_oo
    end

    return trivial_gcd(p, q)

end

function ngcdk(p::P, q::P, k::Int;
               ϵ = 1e-8, ρ = ϵ) where {T, X, P <: PnPolynomial{T,X}}
    m, n = length(p)-1, length(q)-1
    Sₓ = hcat(convmtx(p,k),  convmtx(q, m-n+k))
    F = qr(Sₓ)
    R = UpperTriangular(F.R)
    x = zeros(T, size(Sₓ, 2))
    σ₋₁ = smallest_singular_value!(x, R, ϵ *  sqrt(m - k + 1), norm(Sₓ,2) * ρ)
    @show x
    u, v, w = initial_uvw(Val(:ispossible), k, p, q, x)
    ρₖ, κ = refine_uvw!(u, v, w, p, q, uv, uw)
    return (u=u, v=v, w=w, θ=ρₖ, κ=κ)
end

function trivial_gcd(p::P, q) where {T, X, P <: PnPolynomial{T, X}}
    u, v, w = one(P), p, q
    return (u=u, v=v, w=w, θ=zero(T), κ=typemax(T))
end

# G = [cong(g.c) cong(g.s); -conj(g.s) cong(g.c)]

# A = QR solve by iteration for largest eigenvalue of A^-1
# modifies w in place

# https://arxiv.org/abs/2103.04196
function smallest_singular_value!(w, R::UpperTriangular{T},
                                  atol, # numerical tolerance
                                  rtol# ||A||_oo * eps(T)
                                  ) where {T}

    # Cant' handle singular matrices
    iszero(det(R)) && return zero(T)

    MAXSTEPS = 50

    sⱼ = sⱼ₋₁ = typemax(T)

    rand!(w)
    w ./= norm(w)

    j = 1
    while true
        sⱼ₋₁ = sⱼ
        sⱼ = smallest_singular_value_one_step!(w, R)

        if sⱼ ≤ atol || abs(sⱼ - sⱼ₋₁) ≤ sⱼ * rtol
            break
        end

        j += 1
        j >= MAXSTEPS && break
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
function solve_u(v::P, w, p, q, k) where {T, X, P<:PnPolynomial{T,X}}
    A = [convmtx(v, k+1); convmtx(w, k+1)]
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

    ρ₀, ρ₁ = one(T), residual_error(p,q,uv,uw)

    # storage
    h, β =  u, dot(u,u)  # h = constant * u₀ is used
    A = JF(h, u, v, w)
    Δfβ = Fmp(dot(h, u) - β, p, q, uv, uw)
    Δz = zeros(T, length(u) + length(v) + length(w))

    steps = 0
    @show steps, ρ₁

    minᵢ, Maxᵢ = 3, 20



    ũ, ṽ, w̃ = copy(u), copy(v), copy(w)
    while ρ₁ > 0.0
        steps += 1

        refine_uvw_step!(ũ, ṽ, w̃,
                         Δz, A, Δfβ)

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

    _, R = qr(A) # uggh, need to save. It is computed in refine_uvw_step!
    σ₂ = smallest_singular_value_one_step!(Δz, UpperTriangular(R))
    κ = 1 / σ₂
            @show :use, ρ₁
    return ρ₁, κ

end

# update u,v,w, uv, uw
# computes update step of zₖ₊₁ = zₖ - Δz; Δz = J(zₖ)⁺(fₕ(u,v,w) - [β;p;q])
function refine_uvw_step!(u, v, w,
                          Δz, J⁺, Δfβ)

    qrsolve!(Δz, J⁺, Δfβ)

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
## Δ = h⋅uᵢ - uᵢ⋅uᵢ
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

function qrsolve!(y::Vector{T}, A, b) where {T}
    y .= A \ b
end

# Fast least-squares solver for full column rank Hessenberg-like matrices
# By Andreas Varga
function qrsolve!(y::Vector{T}, A, b) where {T <: Float64}
    # half the time of ldiv!(y, qr(A), b)
    Base.require_one_based_indexing(A)
    m, n = size(A)
    m < n && error("Column dimension exceeds row dimension")
    _, τ = LinearAlgebra.LAPACK.geqrf!(A)
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
    sⱼ = sⱼ₋₁ = typemax(T)

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
