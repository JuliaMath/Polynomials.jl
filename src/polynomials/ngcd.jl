function ngcd(p::P, q::Q,
              args...; kwargs...) where {T,X,P<:AbstractPolynomial{T,X},
                                         S,Y,Q<:AbstractPolynomial{S,Y}}

    degree(p) < 0  && return (u=q,      v=p, w=one(q),  θ=NaN, κ=NaN)
    degree(p) == 0 && return (u=one(q), v=p, w=q,       θ=NaN, κ=NaN)
    degree(q) < 0  && return (u=one(q), v=p, w=zero(q), θ=NaN, κ=NaN)
    degree(q) == 0 && return (u=one(p), v=p, w=q,       θ=NaN, κ=NaN)
    p == q         && return (u=p,v=one(p),  w=one(p),  θ=NaN, κ=NaN)
    Polynomials.assert_same_variable(p,q)
    
    R = promote_type(float(T), float(S))
    ps = R[pᵢ for pᵢ ∈ coeffs(p)]
    qs = R[qᵢ for qᵢ ∈ coeffs(q)]
    chop!(ps); chop!(qs)
    p′ = NGCD.NCPolynomial(ps)
    q′ = NGCD.NCPolynomial(qs)
    if degree(p′) > 5*degree(q′)
        out = NGCD.ngcd′(p′, q′, args...; kwargs...)
    else
        out = NGCD.ngcd(p′, q′; kwargs...)
    end

    𝑷 = promote_type(P,Q)
    u,v,w = convert.(𝑷, (out.u,out.v,out.w))
    (u=u,v=v,w=w, Θ=out.Θ, κ = out.κ)
end


module NGCD
using Polynomials, LinearAlgebra
# non-copying polynomial for performance reasons
struct NCPolynomial{T <: Number, X} <: Polynomials.StandardBasisPolynomial{T, X}
    coeffs::Vector{T}
    function NCPolynomial{T, X}(coeffs::AbstractVector{T}) where {T <: Number, X}
        iszero(coeffs[end]) && throw(ArgumentError("trim your coeffs; $coeffs"))
        new{T,X}(coeffs)
    end
end

NCPolynomial(coeffs::AbstractVector{T}, var=:x) where {T} = NCPolynomial{T,X}(coeffs)
Polynomials.@register NCPolynomial
Base.broadcastable(p::NCPolynomial) = p.coeffs;
Base.ndims(::Type{<:NCPolynomial}) = 1
Base.copyto!(p::NCPolynomial, x) = (copyto!(p.coeffs, x); chop!(p));


     
"""
    ngcd′(p,q)

When degree(p) ≫ degree(q), this uses a early call to `divrem` to bring about commensurate degrees
before calling `ngcd`.
"""
function ngcd′(p::NCPolynomial{T}, q::NCPolynomial{T};
               atol = eps(real(float(T))),
               rtol = atol, 
               satol= atol,
               srtol= rtol,
               kwargs...
               ) where {T}


    a, b = divrem(p,q)

    # check if a=u (p,q) ≈ (aq,q)
    if isapprox(p, a*q, atol=atol, rtol=rtol)
        return (u=a, v=p, w=q, θ=NaN, κ=NaN)
    else
        ngcd(q, b; atol=100atol, rtol=100rtol, kwargs...)
    end
end



function ngcd(p::NCPolynomial{T,X},
              q::NCPolynomial{T,X};
              scale::Bool=false, #norm(ps) > 1e6, #false,
              atol = eps(real(T)),
              rtol = Base.rtoldefault(real(T)),
              satol = eps(real(T))^(5/6),
              srtol = eps(real(T)),
              verbose=false,
              minⱼ = -1
              ) where {T <: AbstractFloat, X}

    m, n = degree.((p, q))
    vw = true
    if m < n
        out = _ngcd(q, p; scale=scale,
                    atol=atol, rtol=rtol, satol=satol,
                    srtol=srtol,
                    verbose=verbose, minⱼ = minⱼ)
        # switch v,w
        return (u=out.u, v=out.w, w=out.v, Θ=out.Θ, κ=out.κ)
    end

    if scale
        p ./= norm(p)
        q ./= norm(q)
    end
    
    # storage
    A0 = zeros(T, m+1, 2)
    A0[:,1] = coeffs(p)
    A0[end-n:end,2] = coeffs(q)

    # pre-allocate storage for Sylvester Matrices, S₁, S₂...
    Q = zeros(T, m + n, m + n)
    R = zeros(T, m + n, m + n)
    Sₓ = hcat(convmtx(p,1),  convmtx(q, m-n+1))
    local x::Vector{T}

    j = n  # We count down Sn, S_{n-1}, ..., S₂, S₁
    
    F = qr(Sₓ)
    nr, nc = size(Sₓ) # m+1, m-n+2
    Q[1:nr, 1:nr] = F.Q
    R[1:nc, 1:nc] = F.R

    while true

        V = view(R, 1:nc, 1:nc)
        flag, σ, x = smallest_singular_value(V, satol *  sqrt(1 + m - j), srtol)
        verbose && println("------ degree $j ----- σ₁: $σ  --- $flag")

        if (flag == :iszero || flag == :ispossible)

            u, v, w = initial_uvw(Val(flag), j, p, q, x)
            flag, ρ₁, σ₂, ρ = refine_uvw!(u,v,w, p, q, atol, rtol)

            verbose && println("   --- Θᵏ: $ρ₁ --- $flag (ρ=$(ρ))")

            if flag == :convergence
                return (u=u, v=v, w=w, Θ=ρ₁, κ=σ₂) # (u,v,w) verified
            end
        end
        
        # reduce possible degree of u and try again with Sⱼ₋₁
        # unless we hit specified minimum, in which case return it
        if j == minⱼ
            u, v, w = initial_uvw(Val(:ispossible), j, p, q, x)
            flag, ρ₁, σ₂, ρ = refine_uvw!(u,v,w, p, q, atol, rtol)
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
    flag, ρ₁, κ, ρ = refine_uvw!(u,v,w, p, q, atol, rtol)
    return (u=u, v=v, w=w, Θ=ρ₁, κ=κ)

end

# fix the degree, k
function ngcd(p::P,
              q::P,
              k::Int;
              kwargs...
              ) where {T <: AbstractFloat,X, P <: AbstractPolynomial{T,X}}

    m, n = degree.((p,q))

    if m < n
        out = ngcd(q, p, k, atol=Inf, rtol=Inf)
        return (u=out.u, v=out.w, w=out.v, Θ=out.Θ, κ=out.κ)
    end

    #    u,v,w = initial_uvw(Val(:iszero), k, ps, qs, nothing)
    Sⱼ = [convmtx(p, n-k+1) convmtx(q, m-k+1)]
    F = qr(Sⱼ)
    flag, σ, x = smallest_singular_value(F.R, eps(T) *  sqrt(1 + m - k), eps(T))
    if flag != :iszero
        w, v = P(x[1:n-k+1]), P(-x[n-k+2:end])
        u = solve_u(v,w,p,q,k)
    else
        u,v,w = initial_uvw(Val(flag), k, p, q, nothing)
    end
    flag, ρ₁, κ, ρ = refine_uvw!(u,v,w, p, q, Inf, Inf)
    return (u=u, v=v, w=w, Θ=ρ₁, κ=κ) 

end

# function ngcd(ps::Vector{T}, qs::Vector{S}, k::Int; kwargs...) where {T, S}
#     ps′,qs′ = promote(float.(ps), float(qs))
#     ngcd(ps′, qs′, k; kwargs...)
# end

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

## --------------------------------------------------
## Refine u,v,w

## Find u₀,v₀,w₀ from right singular vector
function initial_uvw(::Val{:ispossible}, j, p::P, q, x) where {T,X,P<:AbstractPolynomial{T,X}}

    # Sk*[w;-v] = 0, so pick out v,w after applying permuation
    m,n = degree.((p, q))
    vᵢ = vcat(2:m-n+2, m-n+4:2:length(x))
    wᵢ = m-n+3 > length(x) ? [1] : vcat(1, (m-n+3):2:length(x))
    v = P(-x[vᵢ])
    w = P(x[wᵢ])
    # p194 3.9 C_k(v) u = p or Ck(w) u = q; this uses 10.2
    u = solve_u(v,w,p,q,j)
    return u,v,w
end

function initial_uvw(::Val{:iszero}, j, p::P, q, x) where {T,X,P<:AbstractPolynomial{T,X}}
    m,n = degree.((p,q))
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

function initial_uvw(::Val{:constant}, j, p::P, q, x) where {T,X,P<:AbstractPolynomial{T,X}}
    u = one(P)
    w = q
    v = p
    u,v,w
end


# find estimate for σ₂, used in a condition number (κ = 1/σ)
function σ₂(J)
    F = qr(J)
    flag, σ, x = smallest_singular_value(F.R)
    σ
end
    
## attempt to refine u,v,w
## check that [u ⊗ v; u ⊗ w] ≈ [p; q]
function refine_uvw!(u::P, v, w, p, q, atol, rtol) where {T,X,P<:AbstractPolynomial{T,X}}
    
    m, n, l =  degree.((u, v, w))

    uv = u * v
    uw = u * w
    ρ₀, ρ₁ = one(T), residual_error(p,q,uv,uw)

    # storage
    A = zeros(T, JF_size(u, v, w)...)
    b = zeros(T, 1 + length(p) + length(q))     #b = zeros(T, 1 + length(p) + length(q))
    Δf = zeros(T, m + n + l + 3)
    steps = 0

    h, β =  u, norm(u)^2
    minᵢ, Maxᵢ = 5, 20
    κ = NaN

    JF!(A, h, u, v, w)
    Fmp!(b,  dot(h,u) - β, p, q, uv, uw)

    while ρ₁ > 0.0

        # update A, b, then solve A\b
        qrsolve!(Δf, A, b)

        # m + n = degree(p)
        # m + l = degree(q)
        # b has length degree(p)+degree(q) + 3
        Δv, Δw, Δu = P(Δf[1:(n+1)]), P(Δf[(n+1) .+ (1:l+1)]), P(Δf[n+1+l+1 + 1:end])
        
        uv, uw =  (u-Δu) * (v-Δv),  (u-Δu) * (w-Δw)
        ρ₀, ρ′ = ρ₁, residual_error(p, q, uv, uw)

        # don't worry about first few, but aftewards each step must be productive
        # though we can have really bad first steps, which we cap
        if  (steps <= Maxᵢ) && (steps <= minᵢ || ρ′ < 0.95 * ρ₀) && (  ρ′ < 100*ρ₁ )
            ρ₁ = ρ′
            u .-= Δu; v .-= Δv; w .-= Δw
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

## Jacobian F(u,v,w) = [p,p'] is J(u,v,w)
function JF_size(u, v, w)

    m, k, j = degree(u), degree(v), degree(w)
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

    du, dv, dw = degree(u), degree(v), degree(w)    
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
    tot = zero(T)
    for (pᵢ, uvᵢ) in zip(p,uv)
        tot += (pᵢ-uvᵢ)^2
    end
    for (qᵢ, uwᵢ) in zip(q, uw)
        tot += (qᵢ-uwᵢ)^2
    end
    sqrt(tot)
end



## --------------------------------------------------
## utils



"""
    convmtx(v, n::Int)
    convmtx!(M,v,n)
    convmtx_size(v,n)

Convolution matrix.
C = convmtx(v,n) returns the convolution matrix C for a vector v. 
If q is a column vector of length n, then C*q is the same as conv(v,q). 

"""
convmtx!
function convmtx!(C, v::AbstractVector{T}, n::Int) where T

    #   Form C as the Toeplitz matrix 
    #   C = Toeplitz([v; zeros(n-1)],[v[1]; zeros(n-1));  put Toeplitz code inline

    nv = length(v)-1

    @inbounds for j = 1:n
        C[j:j+nv,j] = v
    end

    nothing
  
end

function convmtx!(C, v::AbstractPolynomial, n::Int) where T
    convmtx!(C, coeffs(v), n)
end

convmtx_size(v::AbstractPolynomial, n) = (n + degree(v), n)
function convmtx(v::AbstractPolynomial{T}, n::Int) where {T}
    d = degree(v)
    C = zeros(T, (n + d, n))
    convmtx!(C, coeffs(v), n)
    C
end


# solve for u from [v,w] \ [p,q]
function solve_u(v::P,w,p,q, k) where {T,X,P<:NCPolynomial{T,X}}
    A = [convmtx(v,k+1); convmtx(w, k+1)]
    b = vcat(coeffs(p), coeffs(q))
    u = P(A \ b)
    return u
end

end
