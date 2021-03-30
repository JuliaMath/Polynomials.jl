abstract type StandardBasisPolynomial{T,X} <: AbstractPolynomial{T,X} end

function showterm(io::IO, ::Type{<:StandardBasisPolynomial}, pj::T, var, j, first::Bool, mimetype) where {T}

    if pj === zero(T) return false end

    pj = printsign(io, pj, first, mimetype)

    if !(pj === one(T) && !(showone(T) || j == 0))
        printcoefficient(io, pj, j, mimetype)
    end

    printproductsign(io, pj, j, mimetype)
    printexponent(io, var, j, mimetype)
    return true
end

"""
    evalpoly(x, p::StandardBasisPolynomial)
    p(x)

Evaluate the polynomial using [Horner's Method](https://en.wikipedia.org/wiki/Horner%27s_method), also known as synthetic division, as implemented in `evalpoly` of base `Julia`.

# Examples
```jldoctest
julia> using Polynomials

julia> p = Polynomial([1, 0, 3])
Polynomial(1 + 3*x^2)

julia> p(0)
1

julia> p.(0:3)
4-element Array{Int64,1}:
  1
  4
 13
 28
```
"""
evalpoly(x, p::StandardBasisPolynomial) = EvalPoly.evalpoly(x, p.coeffs) # allows  broadcast  issue #209
constantterm(p::StandardBasisPolynomial) = p[0]

domain(::Type{<:StandardBasisPolynomial}) = Interval(-Inf, Inf)
mapdomain(::Type{<:StandardBasisPolynomial}, x::AbstractArray) = x

function Base.convert(P::Type{<:StandardBasisPolynomial}, q::StandardBasisPolynomial)
    if isa(q, P)
        return q
    else
        T = _eltype(P,q)
        X = indeterminate(P,q)
        return ⟒(P){T,X}([q[i] for i in 0:degree(q)])
    end
end

function Base.one(::Type{P}) where {P<:StandardBasisPolynomial}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}(ones(T,1))
end
function variable(::Type{P}) where {P<:StandardBasisPolynomial}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}([zero(T),one(T)])
end

## multiplication algorithms for computing p * q.
## default multiplication between same type.
## subtypes might relax to match T,S to avoid one conversion
function Base.:*(p::P, q::P) where {T,X, P<:StandardBasisPolynomial{T,X}}
    cs = ⊗(P, coeffs(p), coeffs(q))
    P(cs)
end
                                    
## put here, not with type defintion, in case reuse is possible
function ⊗(P::Type{<:StandardBasisPolynomial}, p::Vector{T}, q::Vector{S}) where {T,S}
    R = promote_type(T,S)
    fastconv(convert(Vector{R}, p), convert(Vector{R},q))
end

## Static size of product makes generated functions  a good choice
## from https://github.com/tkoolen/StaticUnivariatePolynomials.jl/blob/master/src/monomial_basis.jl
## convolution of two tuples
@generated function ⊗(::Type{<:StandardBasisPolynomial}, p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}
    P = M + N - 1
    exprs = Any[nothing for i = 1 : P]
    for i in 1 : N
        for j in 1 : M
            k = i + j - 1
            if exprs[k] === nothing
                exprs[k] = :(p1[$i] * p2[$j])
            else
                exprs[k] = :(muladd(p1[$i], p2[$j], $(exprs[k])))
            end
        end
    end

    return quote
        Base.@_inline_meta
        tuple($(exprs...))        
    end

end

function ⊗(P::Type{<:StandardBasisPolynomial}, p::Dict{Int,T}, q::Dict{Int,S}) where {T,S}
    R = promote_type(T,S)
    c = Dict{Int,R}()
    for (i,pᵢ) ∈ pairs(p)
        for (j,qⱼ) ∈ pairs(q)
            cᵢⱼ = get(c, i+j, zero(R))
            @inbounds c[i+j] = muladd(pᵢ, qⱼ, cᵢⱼ)
        end
    end
    c
end

## ---
function fromroots(P::Type{<:StandardBasisPolynomial}, r::AbstractVector{T}; var::SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j in 1:n, i in j:-1:1
        c[(i + 1)] = c[(i + 1)] - r[j] * c[i]
    end
    #return P(c, var)
    if sort(r[imag(r).>0], lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y)) == sort(conj(r[imag(r).<0]), lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y))
        c = real(c)             # if complex poles come in conjugate pairs, the coeffs are real
    end
    return P(reverse(c), var)
end


function derivative(p::P, order::Integer = 1) where {T, X, P <: StandardBasisPolynomial{T, X}}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))

    # we avoid usage like Base.promote_op(*, T, Int) here, say, as
    # Base.promote_op(*, Rational, Int) is Any, not Rational in analogy to
    # Base.promote_op(*, Complex, Int)
    R = eltype(one(T)*1)
    Q = ⟒(P){R,X}
    
    order == 0 && return p
    hasnan(p) && return Q(R[NaN])
    order > length(p) && return zero(Q)
    d = degree(p)
    d <= 0 && return zero(Q)
    n = d + 1
    a2 = Vector{R}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    Q(a2)
end

function integrate(p::P) where {T, X, P <: StandardBasisPolynomial{T, X}}
    R = eltype(one(T)/1)
    Q = ⟒(P){R,X}    

    hasnan(p) && return Q([NaN])
    iszero(p) && return zero(Q)

    n = length(p)
    as = Vector{R}(undef, n + 1)
    as[1] = zero(R)
    for (i, pᵢ) ∈ pairs(p)
        i′ = i + 1
        @inbounds as[i′+1] = pᵢ/i′
    end
    return Q(as)
end

function Base.divrem(num::P, den::Q) where {T, P <: StandardBasisPolynomial{T}, S, Q <: StandardBasisPolynomial{S}}

    assert_same_variable(num, den) 
    X = indeterminate(num)


    n = degree(num)
    m = degree(den)

    m == -1 && throw(DivideError())
    if m == 0 && den[0] ≈ 0 throw(DivideError()) end

    R = eltype(one(T)/one(S))

    deg = n - m + 1

    if deg ≤ 0
        return zero(P), num
    end

    q_coeff = zeros(R, deg)
    r_coeff = R[ num[i-1] for i in 1:n+1 ]

    @inbounds for i in n:-1:m
        q = r_coeff[i + 1] / den[m]
        q_coeff[i - m + 1] = q
        @inbounds for j in 0:m
            elem = den[j] * q
            r_coeff[i - m + j + 1] -= elem
        end
    end
    resize!(r_coeff, min(length(r_coeff), m))

    return _convert(num, q_coeff), _convert(num, r_coeff)

end

"""
    gcd(p1::StandardBasisPolynomial, p2::StandardBasisPolynomial; method=:eculidean, kwargs...)

Find the greatest common divisor.

By default, uses the Euclidean division algorithm (`method=:euclidean`), which is susceptible to floating point issues.

Passing `method=:noda_sasaki` uses scaling to circumvent some of these.

Passing `method=:numerical` will call the internal method `NGCD.ngcd` for the numerical gcd. See the help page of [`Polynomials.NGCD.ngcd`](@ref) for details.
"""
function Base.gcd(p1::P, p2::Q, args...;
                  method=:euclidean,
                  kwargs...
                  ) where {T, P <: StandardBasisPolynomial{T}, Q <: StandardBasisPolynomial{T}}

    gcd(Val(method), p1, p2, args...; kwargs...)
end

function Base.gcd(::Val{:euclidean},
                  r₀::StandardBasisPolynomial{T}, r₁;
                  atol=zero(real(T)),
                  rtol=Base.rtoldefault(real(T)),
                  kwargs...) where {T}


    iter = 1
    itermax = length(r₁)

    while !iszero(r₁) && iter ≤ itermax
        _, rtemp = divrem(r₀, r₁)
        r₀ = r₁
        r₁ = truncate(rtemp; atol=atol, rtol=rtol)
        iter += 1
    end
    return r₀
end


Base.gcd(::Val{:noda_sasaki}, p, q; kwargs...) = gcd_noda_sasaki(p,q; kwargs...)

"""
     gcd_noda_sasaki(p,q; atol,  rtol)

Greatest common divisor of two polynomials.
Compute the greatest common divisor `d` of two polynomials `a` and `b` using
the Euclidian Algorithm with scaling as of [1].

References:

[1] M.-T. Noda and T. Sasaki. Approximate GCD and its application to ill-conditioned  algebraic equations. J. Comput. Appl. Math., 38:335–351, 1991.

Author: Andreas Varga

Note: requires Julia `v1.2` or greater.
"""
function  gcd_noda_sasaki(p::StandardBasisPolynomial{T}, q::StandardBasisPolynomial{S};
                          atol::Real=zero(real(promote_type(T,S))),
                          rtol::Real=Base.rtoldefault(real(promote_type(T,S)))
                          ) where {T,S}
    ⟒(typeof(p)) == ⟒(typeof(q)) ||  return gcd_noda_sasaki(promote(p,q);  atol=atol, rtol=rtol)
    ## check symbol
    a, b = coeffs(p), coeffs(q)
    as =  _gcd_noda_sasaki(a,b, atol=atol,  rtol=rtol)

    _convert(p, as)
end

function _gcd_noda_sasaki(a::Vector{T}, b::Vector{S};
              atol::Real=zero(real(promote_type(T,S))),
              rtol::Real=Base.rtoldefault(real(promote_type(T,S)))
              ) where {T,S}

    R = eltype(one(T)/one(S))

    na1 = findlast(!iszero,a) # degree(a) + 1
    na1 === nothing && return(ones(R, 1))

    nb1 = findlast(!iszero,b) # degree(b) + 1
    nb1 === nothing && return(ones(R, 1))

    a1 = R[a[i] for i in 1:na1]
    b1 = R[b[i] for i in 1:nb1]
    a1 ./= norm(a1)
    b1 ./= norm(b1)

    tol = atol + rtol

    # determine the degree of GCD as the nullity of the Sylvester matrix
    # this computation can be replaced by simply setting nd = 1, in which case the Sylvester matrix is not formed

    nd = na1 + nb1 - 2 - rank([NGCD.convmtx(a1,nb1-1) NGCD.convmtx(b1,na1-1)], atol = tol) # julia 1.1
    nd == 0 && (return [one(R)])

    sc = one(R)
    while na1 > nd
         na1 < nb1 && ((a1, b1, na1, nb1) = (b1, a1, nb1, na1))
         @inbounds for i in na1:-1:nb1
            s = -a1[i] / b1[nb1]
            sc = max(sc, abs(s))
            @inbounds for j in 1:nb1-1
                a1[i-nb1+j] += b1[j] * s
            end
        end
        a1 ./= sc
        na1 = findlast(t -> (abs(t) > tol),a1[1:nb1-1])
        na1 === nothing && (na1 = 0)
        resize!(a1, na1)
    end

    return nb1 == 1 ? [one(R)] : b1

end


Base.gcd(::Val{:numerical}, p, q, args...; kwargs...) = ngcd(p,q, args...; kwargs...).u


## --------------------------------------------------

function companion(p::P) where {T, P <: StandardBasisPolynomial{T}}
    d = length(p) - 1
    d < 1 && throw(ArgumentError("Series must have degree greater than 1"))
    d == 1 && return diagm(0 => [-p[0] / p[1]])


    R = eltype(one(T)/one(T))

    comp = diagm(-1 => ones(R, d - 1))
    ani = 1 / p[end]
    for j in  0:(degree(p)-1)
        comp[1,(d-j)] = -p[j] * ani # along top row has smaller residual than down column
    end
    return comp
end

function  roots(p::P; kwargs...)  where  {T, P <: StandardBasisPolynomial{T}}

    R = eltype(one(T)/one(T))
    d = degree(p)
    if d < 1
        return R[]
    end
    d == 1 && return R[-p[0] / p[1]]

    as = [p[i] for i in 0:d]
    K  = findlast(!iszero, as)
    if K == nothing
        return R[]
    end
    k =  findfirst(!iszero, as)

    k  == K && return zeros(R, k-1)

    # find eigenvalues of the  companion matrix
    comp  = companion(⟒(P)(as[k:K], indeterminate(p)))
    L = eigvals(comp; kwargs...)
    append!(L, zeros(eltype(L), k-1))

    L
end

function vander(P::Type{<:StandardBasisPolynomial}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    @inbounds for i in 1:n
        A[:, i + 1] = A[:, i] .* x
    end
    return A
end


## as noted at https://github.com/jishnub/PolyFit.jl, using method from SpecialMatrices is faster
## https://github.com/JuliaMatrices/SpecialMatrices.jl/blob/master/src/vandermonde.jl
## This is Algorithm 2 of https://www.maths.manchester.ac.uk/~higham/narep/narep108.pdf
## Solve V(αs)⋅x = y where V is (1+n) × (1+n) Vandermonde matrix (Vᵀ in the paper)
function solve_vander!(ys, αs)
    n = length(ys) - 1
    for k in 0:n-1
        for j in n:-1:k+1
            ys[1+j] = (ys[1+j] - ys[1+j-1])/(αs[1+j] - αs[1+j-k-1])
        end
    end
    for k in n-1:-1:0
        for j in k:n-1
            ys[1+j] = ys[1+j] - αs[1+k] * ys[1 + j + 1]
        end
    end

    nothing
end

# intercept one (typical) case for a faster variant
function fit(P::Type{<:StandardBasisPolynomial},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::Integer = length(x) - 1;
             weights = nothing,
             var = :x,) where {T}

    if deg == length(x) -1 && weights == nothing
        _polynomial_fit(P, x, y; var=var)
    else
        _fit(P, x, y, deg; weights=weights, var=var)
    end
end

function _polynomial_fit(P::Type{<:StandardBasisPolynomial}, x::AbstractVector{T}, y; var=:x) where {T}
    R = float(T)
    coeffs = Vector{R}(undef, length(x))
    copyto!(coeffs, y)
    solve_vander!(coeffs, x)
    P(R.(coeffs), var)
end

## --------------------------------------------------
## More accurate fitting for higher degree polynomials

"""
    polyfitA(x, y, n=length(x)-1; var=:x)
    fit(ArnoldiFit, x, y, n=length(x)-1; var=:x)

Fit a degree ``n`` or less polynomial through the points ``(x_i, y_i)`` using Arnoldi orthogonalization of the Vandermonde matrix.

The use of a Vandermonde matrix to fit a polynomial to data is exponentially ill-conditioned for larger values of ``n``. The Arnoldi orthogonalization fixes this problem.

# Returns

Returns an instance of `ArnoldiFit`. This object can be used to evaluate the polynomial. To manipulate the polynomial, the object can be `convert`ed to other polynomial types, though there may be some loss in accuracy when doing polynomial evaluations afterwards for higher-degree polynomials.

# Citations:

The two main functions are translations from example code in:

*VANDERMONDE WITH ARNOLDI*;
PABLO D. BRUBECK, YUJI NAKATSUKASA, AND LLOYD N. TREFETHEN;
[arXiv:1911.09988](https://people.maths.ox.ac.uk/trefethen/vander_revised.pdf)

# Examples:

```
f(x) = 1/(1 + 25x^2)
N = 80; xs = [cos(j*pi/N) for j in N:-1:0];
p = fit(Polynomial, xs, f.(xs));
q = fit(ArnoldiFit, xs, f.(xs));
maximum(abs, p(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 3.304586010148457e16
maximum(abs, q(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 1.1939520722092922e-7

N = 250; xs = [cos(j*pi/N) for j in N:-1:0];
p = fit(Polynomial, xs, f.(xs));
q = fit(ArnoldiFit, xs, f.(xs));
maximum(abs, p(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 3.55318186254542e92
maximum(abs, q(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 8.881784197001252e-16

p = fit(Polynomial, xs, f.(xs), 10); # least-squares fit
q = fit(ArnoldiFit, xs, f.(xs), 10);
maximum(abs, q(x) - p(x) for x ∈ range(-1,stop=1,length=500)) # 4.6775083806238626e-14
Polynomials.norm(q-p, Inf) # 2.2168933355715126e-12 # promotes `q` to `Polynomial`
```

"""
function polyfitA(x, y, n=length(x)-1; var=:x)
    m = length(x)
    T = eltype(y)
    Q = ones(T, m, n+1)
    #H = UpperHessenberg(zeros(T, n+1, n))
    H = zeros(T, n+1, n)

    q = zeros(T, m)

    @inbounds for k = 1:n
        q .= x .* Q[:,k]
        for j in 1:k
            λ = dot(Q[:,j], q)/m
            H[j,k] = λ
            q .-= λ * Q[:,j]
        end
        H[k+1,k] = norm(q)/sqrt(m)
        Q[:,k+1] .= q/H[k+1,k]
    end
    d = Q \ y
    ArnoldiFit{eltype(d),typeof(H),Symbol(var)}(d, H)
end

function polyvalA(d, H::AbstractMatrix{S}, s::T) where {T, S}
    R = promote_type(T,S)
    n = length(d) - 1
    W = ones(R, n+1)
    @inbounds for k in 1:n
        w = s .* W[k]
        for j in 1:k
            w -= H[j,k] * W[j]
        end
        W[k+1] = w/H[k+1,k]
    end
    sum(W[i]*d[i] for i in eachindex(d))
end

# Polynomial Interface
"""
    ArnoldiFit

A polynomial type produced through fitting a degree ``n`` or less polynomial to data ``(x_1,y_1),…,(x_N, y_N), N ≥ n+1``, This uses Arnoldi orthogonalization to avoid the exponentially ill-conditioned Vandermonde polynomial. See [`Polynomials.polyfitA`](@ref) for details.
"""
struct ArnoldiFit{T, M<:AbstractArray{T,2}, X}  <: AbstractPolynomial{T,X}
    coeffs::Vector{T}
    H::M
end
export ArnoldiFit
@register ArnoldiFit
domain(::Type{<:ArnoldiFit}) = Interval(-Inf, Inf)

Base.show(io::IO, mimetype::MIME"text/plain", p::ArnoldiFit) = print(io, "ArnoldiFit of degree $(length(p.coeffs)-1)")

evalpoly(x, p::ArnoldiFit) = polyvalA(p.coeffs, p.H, x)

fit(::Type{ArnoldiFit}, x::AbstractVector{T}, y::AbstractVector{T}, deg::Int=length(x)-1;  var=:x, kwargs...) where{T} = polyfitA(x, y, deg; var=var)

Base.convert(::Type{P}, p::ArnoldiFit) where {P <: AbstractPolynomial} = p(variable(P,indeterminate(p)))




## --------------------------------------------------
"""
    compensated_horner(p::P, x)
    compensated_horner(ps, x)

Evaluate `p(x)` using a compensation scheme of S. Graillat, Ph. Langlois, N. Louve [Compensated Horner Scheme](https://cadxfem.org/cao/Compensation-horner.pdf). Either a `Polynomial` `p` or its coefficients may be passed in.

The Horner scheme has relative error given by

`|(p(x) - p̂(x))/p(x)| ≤ α(n) ⋅ u ⋅ cond(p, x)`, where `u` is the precision (`2⁻⁵³ = eps()/2`)).

The compensated Horner scheme has relative error bounded by

`|(p(x) - p̂(x))/p(x)| ≤  u +  β(n) ⋅ u² ⋅ cond(p, x)`.

As noted, this reflects the accuracy of a backward stable computation performed in doubled working precision `u²`. (E.g., polynomial evaluation of a `Polynomial{Float64}` object through `compensated_horner` is as accurate as evaluation of a `Polynomial{Double64}` object (using the `DoubleFloat` package), but significantly faster.

Pointed out in https://discourse.julialang.org/t/more-accurate-evalpoly/45932/5.
"""
function compensated_horner(p::P, x) where {T, P <: StandardBasisPolynomial{T}}
    compensated_horner(coeffs(p), x)
end

# rexpressed from paper to compute horner_sum in same pass
# sᵢ -- horner sum
# c -- compensating term
@inline function compensated_horner(ps, x)
    n, T = length(ps), eltype(ps)
    aᵢ = ps[end]
    sᵢ = aᵢ * EvalPoly._one(x)
    c = zero(T) * EvalPoly._one(x)
    for i in n-1:-1:1
	aᵢ = ps[i]
        pᵢ, πᵢ = two_product_fma(sᵢ, x)
	sᵢ, σᵢ = two_sum(pᵢ, aᵢ)
        c = fma(c, x, πᵢ + σᵢ)
    end
    sᵢ + c
end

function compensated_horner(ps::Tuple, x::S) where {S}
    ps == () && return zero(S)
    if @generated
        n = length(ps.parameters)
        sσᵢ =:(ps[end] * EvalPoly._one(x), zero(S))
        c = :(zero(S) * EvalPoly._one(x))
        for i in n-1:-1:1
            pπᵢ = :(two_product_fma($sσᵢ[1], x))
	    sσᵢ = :(two_sum($pπᵢ[1], ps[$i]))
            Σ = :($pπᵢ[2] + $sσᵢ[2])
            c = :(fma($c, x, $Σ))
        end
        s = :($sσᵢ[1] + $c)
        s
    else
        compensated_horner(ps, x)
    end
end

# error-free-transformations (EFT) of a∘b = x + y where x, y floating point, ∘ mathematical
# operator (not ⊕ or ⊗)
@inline function two_sum(a, b)
    x = a + b
    z = x - a
    y = (a - (x-z)) + (b-z)
    x, y
end

# return x, y, floats,  with x + y = a * b
@inline function two_product_fma(a, b)
    x = a * b
    y = fma(a, b, -x)
    x, y
end


# Condition number of a standard basis polynomial
# rule of thumb: p̂ a compute value
# |p(x) - p̃(x)|/|p(x)| ≤ α(n)⋅u ⋅ cond(p,x), where u = finite precision of compuation (2^-p)
function LinearAlgebra.cond(p::P, x) where {P <: StandardBasisPolynomial}
    p̃ = map(abs, p)
    p̃(abs(x))/ abs(p(x))
end
