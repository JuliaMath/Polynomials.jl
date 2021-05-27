"""
    ngcd(p,q, [k]; kwargs...)

Find numerical GCD of polynomials `p` and `q`. Refer to [`NGCD.ngcd`](@ref) for details.


In the case `degree(p) ≫ degree(q)`,  a heuristic is employed to first call one step of the Euclidean gcd approach, and then call `ngcd` with relaxed tolerances.

"""
function ngcd(p::P, q::Q,
              args...; kwargs...) where {T,X,P<:StandardBasisPolynomial{T,X},
                                         S,Y,Q<:StandardBasisPolynomial{S,Y}}

    if (degree(q) > degree(p))
        u,w,v,Θ,κ =  ngcd(q,p,args...;kwargs...)
        return (u=u,v=v,w=w, Θ=Θ, κ=κ)
    end
    if degree(p) > 5*(1+degree(q))
        a,b = divrem(p,q)
        return ngcd(q,b, args...; λ=100, kwargs...)
    end

    # easy cases
    degree(p) < 0  && return (u=q,      v=p, w=one(q),  θ=NaN, κ=NaN)
    degree(p) == 0 && return (u=one(q), v=p, w=q,       θ=NaN, κ=NaN)
    degree(q) < 0  && return (u=one(q), v=p, w=zero(q), θ=NaN, κ=NaN)
    degree(q) == 0 && return (u=one(p), v=p, w=q,       θ=NaN, κ=NaN)
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
    out = NGCD.ngcd(p′, q′, args...; kwargs...)

    𝑷 = Polynomials.constructorof(promote_type(P,Q)){R,X}
    u,v,w = convert.(𝑷, (out.u,out.v,out.w))
    if nz > 1
        u *= variable(u)^(nz-1)
    end
    (u=u,v=v,w=w, Θ=out.Θ, κ = out.κ)
    
end

"""
    square_free(p)

Use `ngcd` to identify the square-free part of the polynomial `p`.
"""
square_free(p) = ngcd(p, derivative(p)).v

## ---- the work is done in this module

module NGCD
using Polynomials, LinearAlgebra
import Polynomials: PnPolynomial, constructorof

"""
    ngcd(ps::PnPolynomial{T,X}, qs::PnPolynomial{T,X}, [k::Int]; scale::Bool=false, atol=eps(T), rtol=eps(T), satol=atol, srtol=rtol)

Return `u, v, w, Θ, κ` where ``u⋅v ≈ ps`` and ``u⋅w ≈ qs`` (polynomial multiplication); `Θ` (`\\Theta[tab]`) is the residual error (``‖ [u⋅v,u⋅w] - [ps,qs] ‖``); and `κ` (`\\kappa[tab]`) is the numerical gcd condition number estimate. When `scale=true`, ``u⋅v ≈ ps/‖ps‖`` and ``u⋅w ≈ qs/‖qs‖``

The numerical GCD problem is defined in [1] (5.4). Let ``(p,q)`` be a
polynomial pair with degree m,n. Let Ρmn be set of all such pairs. Any
given pair of polynomials has an exact greatest common divisor, ``u``, of
degree ``k``, defined up to constant factors. Let ``Ρᵏmn`` be the manifold of
all such ``(p,q)`` pairs with exact gcd of degree ``k``. A given pair ``(p,q)`` with exact gcd of degree ``j`` will
have some distance ``Θᵏ`` from ``Pᵏ``.  For a given threshold ``ϵ > 0`` a numerical GCD
of ``(p,q)`` within ``ϵ`` is an exact GCD of a pair ``(p̂,q̂)`` in ``Ρᵏ`` with 

``‖ (p,q) - (p̂,q̂) ‖ <= Θᵏ``, where ``k`` is the largest value for which ``Θᵏ < ϵ``. 

(In the ``ϵ → 0`` limit, this would be the exact GCD.)


Suppose ``(p,q)`` is an ``ϵ`` pertubation from ``(p̂,q̂)`` where ``(p̂,q̂)`` has an exact gcd of degree ``k``, then ``degree(gcdₑ(p,q)) = k``; as ``ϵ → 0``, ``gcdₑ(p,q) → gcd(p̂, q̂)``, and

``limsup_{(p,q)→(p̂,q̂)} inf{ ‖ (u,v,w) - (û,v̂,ŵ) ‖} / ‖ (p,q) - (p̂,q̂) ‖ < κₑ(p,q)``.

``κ`` is called the numerical GCD condition number.


The Zeng algorithm proposes a degree for ``u`` and *if* a triple ``(u,v,w)`` with ``u`` of degree ``k`` and ``(u⋅v, u⋅w)`` in ``Ρᵏmn`` can be found satisfying ``‖ (u⋅v, u⋅w) - (p,q) ‖ < ϵ`` then ``(u,v,w)`` is returned; otherwise the proposed degree is reduced and the process repeats. If not terminated, at degree ``0`` a constant gcd is returned.

The initial proposed degree is the first ``j``,  ``j=min(m,n):-1:1``, where ``Sⱼ`` is believed to have a singular value of ``0`` (``Sⱼ`` being related to the Sylvester matrix of `ps` and `qs`). The verification of the proposed degree is done using a Gauss-Newton iteration scheme holding the degree of ``u`` constant. 

## Scaling:

If `scale=true` the gcd of ``p/‖p‖`` and ``q/‖q‖`` is identified. Scaling can reduce the condition numbers significantly.

## Tolerances:

There are two places where tolerances are utilized:

* in the identification of the rank of ``Sⱼ`` a value ``σ₁ = ‖Rx‖`` is identified. To test if this is zero a tolerance of `max(satol, ‖R‖ₒₒ ⋅ srtol)` is used.

* to test if ``(u ⋅ v, u ⋅ w) ≈ (p,q)`` a tolerance of `max(atol, λ⋅rtol)` is used with `λ` chosen to be  ``(‖(p,q)‖⋅n⋅m)⋅κ′⋅‖A‖ₒₒ`` to reflect the scale of ``p`` and ``q`` and an estimate for the condition number of ``A`` (a Jacobian involved in the solution). 


This seems to work well for a reasaonable range of polynomials, however there can be issues: when the degree of ``p`` is much larger than the degree of ``q``, these choices can fail; when a higher rank is proposed, then too large a tolerance for `rtol` or `atol` can lead to a false verification; when a tolerance for `atol` or `rtol` is too strict, then a degree may not be verified. 

!!! note:
    These tolerances are adjusted from those proposed in [1].

## Specified degree:

When `k` is specified, a value for ``(u,v,w)`` is identified with ``degree(u)=k``. No tolerances are utilized in computing ``Θᵏ``.



Output:

The function outputs a named tuple with names (`u`, `v`, `w`, `Θ`, `κ`). The components `u`,`v`,`w` estimate the gcd and give the divisors. The value `Θ` estimates ``Θᵏ`` and `κ` estimates the numerical condition number.

Example:

```
using Polynomials
x = variable(Polynomial{Float64})
p = (x+10)*(x^9 + x^8/3 + 1)
q = (x+10)*(x^9 + x^8/7 - 6/7)
gcd(p,q) # u a constant
gcd(p,q, method=:numerical)  # u a degree 1 polynomial
Polynomials.NGCD.ngcd(coeffs(p), coeffs(q), verbose=true) # to see some computations
```

Reference:

[1] The Numerical Greatest Common Divisor of Univariate Polynomials
by Zhonggang Zeng;
[url](http://homepages.neiu.edu/~zzeng/uvgcd.pdf);
[doi](https://doi.org/10.1090/conm/556/11014)

Note: Based on work by Andreas Varga; Requires `VERSION >= v"1.2"`.

"""
function ngcd(p::PnPolynomial{T,X},
              q::PnPolynomial{T,X};
              scale::Bool=false, 
              atol = eps(real(T)),
              rtol = Base.rtoldefault(real(T)),
              satol = eps(real(T))^(5/6),
              srtol = eps(real(T)),
              verbose=false,
              minⱼ = -1,
              λ = 1
              ) where {T, X}

    m,n = length(p)-1, length(q)-1
    @assert m >= n

    ## --- begin
    if scale
        p ./= norm(p)
        q ./= norm(q)
    end
    atol *= λ
    rtol *= λ

    # storage
    A0 = zeros(T, m+1, 2)
    A0[:,1] = coeffs(p)
    A0[end-n:end,2] = coeffs(q)

    # pre-allocate storage for Sylvester Matrices, S₁, S₂...
    Q = zeros(T, m + n, m + n)
    R = zeros(T, m + n, m + n)
    Sₓ = hcat(convmtx(p,1),  convmtx(q, m-n+1))

    uv = copy(p)
    uw = copy(q)
    
    local x::Vector{T}

    F = qr(Sₓ) 
    nr, nc = size(Sₓ) # m+1, m-n+2
    Q[1:nr, 1:nr] .= F.Q
    R[1:nc, 1:nc] .= F.R

    j = n  # We count down Sn, S_{n-1}, ..., S₂, S₁

    while true

        V = view(R, 1:nc, 1:nc)
        flag, σ, x = smallest_singular_value(V, satol *  sqrt(1 + m - j), srtol)
        verbose && println("------ degree $j ----- σ₁: $σ  --- $flag")

        if (flag == :iszero || flag == :ispossible)
            u, v, w = initial_uvw(Val(flag), j, p, q, x)
            flag, ρ₁, σ₂, ρ = refine_uvw!(u,v,w, p, q, uv, uw, atol, rtol)

            verbose && println("   --- Θᵏ: $ρ₁ --- $flag (ρ=$(ρ))")

            if flag == :convergence
                return (u=u, v=v, w=w, Θ=ρ₁, κ=σ₂) # (u,v,w) verified
            end
        end
        
        # reduce possible degree of u and try again with Sⱼ₋₁
        # unless we hit specified minimum, in which case return it
        if j == minⱼ
            u, v, w = initial_uvw(Val(:ispossible), j, p, q, x)
            flag, ρ₁, σ₂, ρ = refine_uvw!(u,v,w, p, q, uv, uw, atol, rtol)
            return (u=u, v=v, w=w, Θ=ρ₁, κ=σ₂)
        end

        j -= 1
        nr += 1
        nc += 2
        nc > nr && break
        extend_QR!(Q,R, nr, nc, A0) # before Q⋅R = Sⱼ, now Q⋅R = Sⱼ₋₁


    end

    # u is a constant
    verbose && println("------ GCD is constant ------")

    u, v, w = initial_uvw(Val(:constant), j, p, q, x)
    flag, ρ₁, κ, ρ = refine_uvw!(u,v,w, p, q, uv, uw, atol, rtol)
    return (u=u, v=v, w=w, Θ=ρ₁, κ=κ)

end

# fix the degree, k
function ngcd(p::P,
              q::P,
              k::Int;
              kwargs...
              ) where {T <: AbstractFloat,X, P <: PnPolynomial{T,X}}

    m::Int, n::Int = length(p)-1, length(q)-1

    #    u,v,w = initial_uvw(Val(:iszero), k, ps, qs, nothing)
    Sⱼ = [convmtx(p, n-k+1) convmtx(q, m-k+1)]
    F = qr(Sⱼ)
    flag, σ, x = smallest_singular_value(F.R, eps(T) *  sqrt(1 + m - k), eps(T))
    u,v,w = initial_uvw(Val(:k), flag, k, p, q, x)
    flag, ρ₁, κ, ρ = refine_uvw!(u,v,w, copy(p), copy(q), copy(p), copy(q),
                                 T(Inf), T(Inf)) # no tolerances
    return (u=u, v=v, w=w, Θ=ρ₁, κ=κ) 

end


## -----

# return guess at smallest singular value and right sinuglar value, x
# for an upper triangular matrix, V
function smallest_singular_value(V::AbstractArray{T,2},
                                 atol=eps(real(T)),
                                 rtol=zero(real(T))) where {T}

    R = UpperTriangular(V)
    k = size(R)[1]/2
    if iszero(det(R)) 
        return (:iszero, zero(T), T[])
    end

    m,n = size(R)

    # we are testing if ‖Ax‖ ≈ 0
    # If x is a perfect 0, but x is approximated by x' due to round off
    # then ‖A(x-x')‖ <= ‖A‖⋅‖x - x'‖ so we use ‖A‖ as scale factor
    δ = max(atol,  norm(R,Inf) * rtol)

    x = ones(T, n)
    y = zeros(T, m)
    σ₀ = σ₁ = Inf*one(real(T))
    steps, minᵢ = 1, 5
    
    while true
        y .= R' \ x # use iteration to converge on singular value
        x .= R  \ y
        x ./= norm(x,2)
        σ₁ = norm(R * x, 2)

        if (steps <= 50) && (steps <= minᵢ || σ₁ < 0.05*σ₀) # decreasing, keep going
            σ₀ = σ₁
        else
            break
        end
        steps += 1
    end

    if σ₁ < δ
        return (:ispossible, σ₁, x)
    else
        return (:constant, σ₁, x)
    end
    
end


## --------------------------------------------------
## Refine u,v,w

## Find u₀,v₀,w₀ from right singular vector
function initial_uvw(::Val{:ispossible}, j, p::P, q::Q, x) where {T,X,
                                                              P<:PnPolynomial{T,X},
                                                              Q<:PnPolynomial{T,X}}

    # Sk*[w;-v] = 0, so pick out v,w after applying permutation
    m,n = length(p)-1, length(q)-1
    vᵢ = vcat(2:m-n+2, m-n+4:2:length(x))
    wᵢ = m-n+3 > length(x) ? [1] : vcat(1, (m-n+3):2:length(x))
    #    v = 𝑷{m-j}(-x[vᵢ])
    v = P(-x[vᵢ])
    w = P(x[wᵢ])
    # p194 3.9 C_k(v) u = p or Ck(w) u = q; this uses 10.2
    u = solve_u(v,w,p,q,j)
    return u,v,w
    
end

function initial_uvw(::Val{:iszero}, j, p::P, q::Q, x) where {T,X,
                                                              P<:PnPolynomial{T,X},
                                                              Q<:PnPolynomial{T,X}}
    
    m,n = length(p)-1, length(q)-1
    S = [convmtx(p, n-j+1) convmtx(q, m-j+1)]

    F = qr(S)
    R = UpperTriangular(F.R)

    if iszero(det(R))
        x = eigvals(R)[:,1]
    else
        x = ones(T, size(R,2))
        x .= R \ (R' \ (x/norm(x)))
        x ./= norm(x)
    end

    w = P(x[1:n-j+1])
    v = P(-x[(n-j+2):end])

    u = solve_u(v,w,p,q,j)
    return u,v,w
end

function initial_uvw(::Val{:constant}, j, p::P, q, x) where {T,X,P<:PnPolynomial{T,X}}
    u = one(P)
    w = q
    v = p
    u,v,w
end

function initial_uvw(::Val{:k}, flag, k, p::P, q, x) where {T,X,P<:PnPolynomial{T,X}}
    flag == :iszero && return initial_uvw(Val(flag), k, p, q, nothing)
    n = length(q)-1
    w, v = P(x[1:n-k+1]), P(-x[n-k+2:end])
    u = solve_u(v,w,p,q,k)
    return (u,v,w)
end
    


# find estimate for σ₂, used in a condition number (κ = 1/σ)
function σ₂(J)
    F = qr(J)
    flag, σ, x = smallest_singular_value(F.R)
    σ
end
    
## attempt to refine u,v,w
## check that [u * v; u * w] ≈ [p; q]
function refine_uvw!(u::U, v::V, w::W, p, q, uv, uw, atol, rtol) where {T,X,
                                                                        U<:PnPolynomial{T,X},
                                                                        V<:PnPolynomial{T,X},
                                                                        W<:PnPolynomial{T,X}}
    
    m,n,l = length(u)-1, length(v)-1, length(w)-1

    mul!(uv, u, v)
    mul!(uw, u, w)

    ρ₀, ρ₁ = one(T), residual_error(p,q,uv,uw)

    # storage
    b = zeros(T, (m+n) + (m+l) + 3) # degree(p) + degree(q) + 3 = 1 + length(p) + length(q))
    Δf = zeros(T, m + n + l + 3)
    steps = 0

    h, β =  u, norm(u)^2
    minᵢ, Maxᵢ = 5, 20
    κ = NaN
    A=zeros(T, JF_size(u, v, w)...)
    JF!(A, h, u, v, w)
    Fmp!(b,  dot(h,u) - β, p, q, uv, uw)

    Δvᵢ = 1:(n+1)
    Δwᵢ = (n+1+1):(n+1+l+1)
    Δuᵢ = (n+1+l+1+1):length(Δf)

    while ρ₁ > 0.0

        # update A, b, then solve A\b
        qrsolve!(Δf, A, b)

        # m + n = degree(p)
        # m + l = degree(q)
        # b has length degree(p)+degree(q) + 3
        Δv = view(Δf, Δvᵢ) 
        Δw = view(Δf, Δwᵢ) 
        Δu = view(Δf, Δuᵢ)
        
        u .-= Δu
        v .-= Δv
        w .-= Δw

        mul!(uv, u, v)
        mul!(uv, u, w)
        
        ρ₀, ρ′ = ρ₁, residual_error(p, q, uv, uw)

        # don't worry about first few, but aftewards each step must be productive
        # though we can have really bad first steps, which we cap
        if  (steps <= Maxᵢ) && (steps <= minᵢ || ρ′ < 0.95 * ρ₀) && (  ρ′ < 100*ρ₁ )
            ρ₁ = ρ′
            steps += 1
        else
            break
        end

        # update A,b for next iteration
        JF!(A, h, u, v, w)
        Fmp!(b,  dot(h,u) - β, p, q, uv, uw)
        
    end


    # this is a heuristic
    # sensitivity is Δu / Δp <= ‖ A+ ‖ = κ
    # we use an estimate for ‖(p,q)‖ error times ‖A⁺‖⋅‖A‖ₒₒ
    κ = 1/σ₂(A) # ≈ ‖A⁺‖
    λ = norm((norm(p), norm(q))) * (m * n) * min(1, κ) * norm(A, Inf)
    ρ = max(atol, rtol * λ)

    if ρ₁ <= ρ
        return :convergence, ρ₁, κ, ρ
    else
        return :no_convergence, ρ₁, κ, ρ
    end
    
end

## ---- QR factorization

function qrsolve!(y::Vector{T}, A, b) where {T}
    y .= qr(A) \ b
end

# # Fast least-squares solver for full column rank Hessenberg-like matrices
# # By Andreas Varga
function qrsolve!(y::Vector{T}, A, b) where {T <: Float64}
    Base.require_one_based_indexing(A)
    m, n = size(A) 
    m < n && error("Column dimension exceeds row dimension") 
    _, τ = LinearAlgebra.LAPACK.geqrf!(A)
    T <: Complex ? tran = 'C' : tran = 'T'
    LinearAlgebra.LAPACK.ormqr!('L', tran, A, τ, view(b,:,1:1))
    y .= UpperTriangular(triu(A[1:n,:]))\b[1:n]
end

# extend QR to next size
# Q gets a 1 in nc,nc, 0s should be elswhere
function extend_QR!(Q,R, nr, nc, A0)


    #old Q is m x m, old R is n x n; we add to these
    n = nc-2 
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


## Jacobian F(u,v,w) = [p,p'] is J(u,v,w)
function JF_size(u, v, w)

    m, k, j = length(u)-1, length(v)-1, length(w)-1
    n, l = m + k, m + j

    ai, aj = convmtx_size(v, m + 1)
    bi, bj = convmtx_size(u, k + 1)
    di, dj = convmtx_size(w, m + 1)
    fi, fj = convmtx_size(u, j + 1)
    ci, cj = ai, fj
    ei, ej = di, bj

    (1 + ai + di, aj + bj + cj)
end

# Jacobian matrix
function JF(u::Vector{U}, v::Vector{V}, w::Vector{W}) where {U,V,W}
    R = promote_type(U,V, W)
    M = zeros(R, JF_size(u, v, w)...)
    JF!(M, u, v, w)
    M
end

function JF!(M, h,  u::P, v, w) where {T,X,P<:AbstractPolynomial{T,X}}

    du, dv, dw = length(u)-1, length(v)-1, length(w)-1
    m, n = du + dv, du + dw

    # JF_size should return these
    r11,c11 = convmtx_size(u, dv+1)
    r13,c13 = convmtx_size(v, du+1)
    r22,c22 = convmtx_size(u, dw+1)
    r23,c23 = convmtx_size(w, du+1)

    J11 = view(M, 1:r11, 1:c11)
    J13 = view(M, 1:r13, c11 + c22 .+ (1:c23))
    J22 = view(M, r11 .+ (1:r22), c11 .+ (1:c22))
    J23 = view(M, r13 .+ (1:r23), (c11 + c22) .+ (1:c23))
    convmtx!(J11, u, dv+1)
    convmtx!(J13, v, du+1)
    convmtx!(J22, u, dw+1)
    convmtx!(J23, w, du+1)
    M[end, end-du:end] = coeffs(h)'
    
    return nothing
end

## compute F(u,v,w) - [p, p'] = [u*v, u*w] - [p, p']
function Fmp!(b, γ, p, q, uv, uw)
    b[end] = γ
    for i in 1:1+length(p)-1
        j = i
        b[i] = uv[j] - p[j]
    end
    for i in 1+length(p):length(b)-1
        j = i - length(p)
        b[i] = uw[j] - q[j]
    end
    return nothing
end


function residual_error(p::P,q,uv,uw) where {T,X,P<:AbstractPolynomial{T,X}}
    tot = zero(real(T))
    for (pᵢ, uvᵢ) in zip(p,uv)
        tot += norm(pᵢ-uvᵢ)^2
    end
    for (qᵢ, uwᵢ) in zip(q, uw)
        tot += norm(qᵢ-uwᵢ)^2
    end
    sqrt(tot)
end



## ---- utils

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


# solve for u from [v,w] \ [p,q]
function solve_u(v::P,w,p,q, k) where {T,X,P<:PnPolynomial{T,X}}
    A = [convmtx(v,k+1); convmtx(w, k+1)]
    b = vcat(coeffs(p), coeffs(q))
    u = P(A \ b)
    return u
end

end


