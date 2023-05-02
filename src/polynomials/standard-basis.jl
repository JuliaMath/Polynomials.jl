abstract type StandardBasisPolynomial{T,X} <: AbstractPolynomial{T,X} end
abstract type LaurentBasisPolynomial{T,X} <: StandardBasisPolynomial{T,X} end

function showterm(io::IO, ::Type{<:StandardBasisPolynomial}, pj::T, var, j, first::Bool, mimetype) where {T}

    if _iszero(pj) return false end

    pj = printsign(io, pj, first, mimetype)

    if hasone(T)
        if !(_isone(pj) && !(showone(T) || j == 0))
            printcoefficient(io, pj, j, mimetype)
        end
    else
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
4-element Vector{Int64}:
  1
  4
 13
 28
```
"""
function evalpoly(x, p::StandardBasisPolynomial{T}) where {T}
    # the zero polynomial is a *special case*
    iszero(p) && return zero(x) * zero(T)
    EvalPoly.evalpoly(x, p.coeffs) # allows  broadcast  issue #209
end

constantterm(p::StandardBasisPolynomial) = p[0]

domain(::Type{<:StandardBasisPolynomial}) = Interval{Open,Open}(-Inf, Inf)
mapdomain(::Type{<:StandardBasisPolynomial}, x::AbstractArray) = x

function Base.convert(P::Type{<:StandardBasisPolynomial}, q::StandardBasisPolynomial)
    if isa(q, P)
        return q
    else
        minimumexponent(P) <= minimumexponent(q) ||
            throw(ArgumentError("a $P can not have a minimum exponent of $(minimumexponent(q))"))
        T = _eltype(P,q)
        X = indeterminate(P,q)
        return ⟒(P){T,X}([q[i] for i in eachindex(q)])
    end
end

# treat p as a *vector* of coefficients
Base.similar(p::StandardBasisPolynomial, args...) = similar(coeffs(p), args...)

function Base.one(::Type{P}) where {P<:StandardBasisPolynomial}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}(ones(T,1))
end
function variable(::Type{P}) where {P<:StandardBasisPolynomial}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}([zero(T), one(T)])
end

## ---- arithmetic

# can bypass c*one(P)
Base.:+(p::P, c::T) where {T, X, P<:StandardBasisPolynomial{T, X}} = p + ⟒(P)([c], X)

## multiplication algorithms for computing p * q.
## default multiplication between same type.
## subtypes might relax to match T,S to avoid one conversion
function Base.:*(p::P, q::P) where {T,X, P<:StandardBasisPolynomial{T,X}}
    cs = ⊗(P, coeffs(p), coeffs(q))
    P(cs)
end

function ⊗(P::Type{<:StandardBasisPolynomial}, p::Vector{T}, q::Vector{S}) where {T<:Number,S<:Number}
    R = promote_type(T,S)
    fastconv(convert(Vector{R}, p), convert(Vector{R},q))
end

## put here, not with type definition, in case reuse is possible
## `conv` can be used with matrix entries, unlike `fastconv`
function conv(p::Vector{T}, q::Vector{S}) where {T,S}
    (isempty(p) || isempty(q)) && return promote_type(T, S)[]
    as = [p[1]*q[1]]
    z = zero(eltype(as)) * as[1]
    n,m = length(p)-1, length(q)-1
    for i ∈ 1:n+m
        Σ = z
        for j ∈ max(0, i-m):min(i,n)
            Σ = muladd(p[1+j], q[1 + i-j], Σ)
        end
        push!(as, Σ)
    end
    as
end

function ⊗(P::Type{<:StandardBasisPolynomial}, p::Vector{T}, q::Vector{S}) where {T,S}
    conv(p, q)
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
            if isnothing(exprs[k])
                exprs[k] = :(p1[$i] * p2[$j])
            else
                exprs[k] = :(muladd(p1[$i], p2[$j], $(exprs[k])))
            end
        end
    end

    return quote
        Base.@_inline_meta # 1.8 deprecation
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

## Define a dot product on vectors of polynomials, where polynomials are treated as
## as the scalar types -- the dot is the sum of pairwise products. Might change in
## the future if this causes confusion with the "vector" interpretation of polys.
LinearAlgebra.dot(xv::AbstractArray{T}, yv::AbstractArray{T}) where {T <: StandardBasisPolynomial} =
    sum(conj(x)*y for (x,y) = zip(xv, yv))

## ---

function fromroots(P::Type{<:StandardBasisPolynomial}, r::AbstractVector{T}; var::SymbolLike = Var(:x)) where {T <: Number}
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

## ----

function derivative(p::P, order::Integer = 1) where {T, X, P <: StandardBasisPolynomial{T, X}}

    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return p
    d = degree(p)
    order > d  && return 0*p
    hasnan(p) && return  ⟒(P)(zero(T)/zero(T), X) # NaN{T}

    n = d + 1
    dp = [reduce(*, (i - order + 1):i, init = p[i]) for i ∈ order:d]
    return ⟒(P)(dp, X)

end

function integrate(p::P) where {T, X, P <: StandardBasisPolynomial{T, X}}

    hasnan(p) && return ⟒(P)(NaN, X)
    iszero(p) && return zero(p)/1

    as = [pᵢ/(i+1) for (i, pᵢ) ∈ pairs(p)]
    pushfirst!(as, zero(constantterm(p)))
    return ⟒(P)(as, X)
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

Passing `method=:numerical` will call the internal method `NGCD.ngcd` for the numerical gcd. See the help page of [`Polynomials.NGCD.ngcd(p,q)`](@ref) for details.
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
the Euclidean Algorithm with scaling as of [1].

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
    isnothing(na1) && return(ones(R, 1))

    nb1 = findlast(!iszero,b) # degree(b) + 1
    isnothing(nb1) && return(ones(R, 1))

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
        isnothing(na1) && (na1 = 0)
        resize!(a1, na1)
    end

    return nb1 == 1 ? [one(R)] : b1

end


Base.gcd(::Val{:numerical}, p, q, args...; kwargs...) = ngcd(p,q, args...; kwargs...).u

uvw(p::P, q::Q, args...;
    method=:euclidean,
    kwargs...
    ) where {T, P <: StandardBasisPolynomial{T}, Q <: StandardBasisPolynomial{T}} =
        uvw(Val(method), p, q; kwargs...)

function uvw(::Val{:numerical}, p::P, q::P; kwargs...) where {P <: StandardBasisPolynomial}
    u,v,w,Θ,κ = ngcd(p,q; kwargs...)
    u,v,w
end

function uvw(V::Val{:euclidean}, p::P, q::P; kwargs...) where {P <: StandardBasisPolynomial}
    u = gcd(V,p,q; kwargs...)
    u, p÷u, q÷u
end

function uvw(::Any, p::P, q::P; kwargs...) where {P <: StandardBasisPolynomial}
    throw(ArgumentError("not defined"))
end

# Some things lifted from
# from https://github.com/jmichel7/LaurentPolynomials.jl/blob/main/src/LaurentPolynomials.jl.
# This follows mostly that of intfuncs.jl which is specialized for integers.
function Base.gcdx(a::P, b::P) where {T,X,P<:StandardBasisPolynomial{T,X}}
    # a0, b0=a, b
    s0, s1=one(a), zero(a)
    t0, t1 = s1, s0
    # The loop invariant is: s0*a0 + t0*b0 == a
    x, y = a, b
    while y != 0
        q, r =divrem(x, y)
        x, y = y, r
        s0, s1 = s1, s0 - q*s1
        t0, t1 = t1, t0 - q*t1
    end
    (x, s0, t0)./x[end]
end

Base.gcdx(a::StandardBasisPolynomial{T,X}, b::StandardBasisPolynomial{S,X}) where {S,T,X} =
    gcdx(promote(a, b)...)

# p^m mod q
function Base.powermod(p::StandardBasisPolynomial, x::Integer, q::StandardBasisPolynomial)
    x==0 && return one(q)
    b=p%q
    t=prevpow(2, x)
    r=one(p)
    while true
        if x>=t
            r=(r*b)%q
            x-=t
        end
        t >>>= 1
        t<=0 && break
        r=(r*r)%q
    end
    r
end


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

# block companion matrix
function companion(p::P) where {T, M <: AbstractMatrix{T}, P<:StandardBasisPolynomial{M}}
    C₀, C₁ = companion_pencil(p)
    C₀ * inv(C₁) # could be more efficient
end

# the companion pencil is `C₀`, `C₁` where `λC₁ - C₀` is singular for
# eigen values of the companion matrix: `vᵀ(λC₁ - C₀) = 0` => `vᵀλ = vᵀ(C₀C₁⁻¹)`, where `C₀C₁⁻¹`
# is the companion matrix.
function companion_pencil(p::P) where {T, P<:StandardBasisPolynomial{T}}
    n = degree(p)
    C₁ = diagm(0 => ones(T, n))
    C₁[end,end] = p[end]

    C₀ = diagm(-1 => ones(T, n-1))
    for i ∈ 0:n-1
        C₀[1+i,end] = -p[i]
    end
    C₀, C₁
end

# block version
function companion_pencil(p::P) where {T, M <: AbstractMatrix{T}, P<:StandardBasisPolynomial{M}}

    m,m′ = size(p[0])
    @assert m == m′  # square matrix

    n = degree(p)

    C₀ = zeros(T, n*m, n*m)
    C₁ = zeros(T, n*m, n*m)


    Δ = 1:m
    for i ∈ 1:n-1
        C₁[(i-1)*m .+ Δ, (i-1)*m .+ Δ] .= I(m)
        C₀[i*m .+ Δ, (i-1)*m .+ Δ] .= I(m)
    end
    C₁[(n-1)*m .+ Δ, (n-1)*m .+ Δ] .= p[end]
    for i ∈ 0:n-1
        C₀[i*m .+ Δ, (n-1)*m .+ Δ] .= -p[i]
    end

    C₀, C₁
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
    if isnothing(K)
        return R[]
    end
    k =  findfirst(!iszero, as)

    k  == K && return zeros(R, k-1)

    # find eigenvalues of the companion matrix of the 0-deflated polynomial
    comp  = companion(⟒(P)(as[k:K], indeterminate(p)))
    L = eigvals(comp; kwargs...)
    append!(L, zeros(eltype(L), k-1))

    L
end

function vander(P::Type{<:StandardBasisPolynomial}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    vander(P, x, 0:n)
end

# skip some degrees
function vander(P::Type{<:StandardBasisPolynomial}, x::AbstractVector{T}, degs) where {T <: Number}
    A = Matrix{T}(undef, length(x),  length(degs))
    Aᵢ = one.(x)

    i′ = 1
    for i ∈ 0:maximum(degs)
        if i ∈ degs
            A[:, i′] = Aᵢ
            i′ += 1
        end
        for (i, xᵢ) ∈ enumerate(x)
            Aᵢ[i] *= xᵢ
        end
    end
    A
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

    if deg == length(x) -1 && isnothing(weights)
        _polynomial_fit(P, x, y; var=var)
    else
        _fit(P, x, y, deg; weights=weights, var=var)
    end
end

"""
    fit(P::Type{<:StandardBasisPolynomial}, x, y, J, [cs::Dict{Int, T}]; weights, var)

Using constrained least squares, fit a polynomial of the type
`p = ∑_{i ∈ J} aᵢ xⁱ + ∑ cⱼxʲ` where `cⱼ` are fixed non-zero constants

* `J`: a collection of degrees to find coefficients for
* `cs`: If given, a `Dict` of key/values, `i => cᵢ`, which indicate the degree and value of the fixed non-zero constants.

The degrees in `cs` and those in `J` should not intersect.

Example
```
x = range(0, pi/2, 10)
y = sin.(x)
P = Polynomial
p0 = fit(P, x, y, 5)
p1 = fit(P, x, y, 1:2:5)
p2 = fit(P, x, y, 3:2:5, Dict(1 => 1))
[norm(p.(x) - y) for p ∈ (p0, p1, p2)] # 1.7e-5, 0.00016, 0.000248
```

"""
function fit(P::Type{<:StandardBasisPolynomial},
             x::AbstractVector{T},
             y::AbstractVector{T},
             J,
             cs=nothing;
             weights = nothing,
             var = :x,) where {T}
    _fit(P, x, y, J; weights=weights, var=var)
end


function fit(P::Type{<:StandardBasisPolynomial},
             x::AbstractVector{T},
             y::AbstractVector{T},
             J,
             cs::Dict{Int, S};
             weights = nothing,
             var = :x,) where {T,S}

    for i ∈ J
        haskey(cs, i) && throw(ArgumentError("cs can't overlap with deg"))
    end

    # we subtract off `∑cᵢ xⁱ`ⱼ from `y`;
    # fit as those all degrees not in J are 0,
    # then add back the constant coefficients

    q = SparsePolynomial(cs)
    y′ = y - q.(x)

    p = fit(P, x, y′, J; weights=weights, var=var)

    return p + q
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

This representation is useful for *evaluating* the polynomial, but does not lend itself to other polynomial manipulations.

# Returns

Returns an instance of `ArnoldiFit`. This object can be used to evaluate the polynomial. To manipulate the polynomial, the object can be `convert`ed to other polynomial types, though there may be some loss in accuracy when doing polynomial evaluations afterwards for higher-degree polynomials.

# Citations:

The two main functions are translations from example code in:

*VANDERMONDE WITH ARNOLDI*;
PABLO D. BRUBECK, YUJI NAKATSUKASA, AND LLOYD N. TREFETHEN;
[arXiv:1911.09988](https://people.maths.ox.ac.uk/trefethen/vander_revised.pdf)

For more details, see also:

Lei-Hong Zhang, Yangfeng Su, Ren-Cang Li. Accurate polynomial fitting and evaluation via Arnoldi. Numerical Algebra, Control and Optimization. doi: 10.3934/naco.2023002



# Examples:

```
f(x) = 1/(1 + 25x^2)
N = 80; xs = [cos(j*pi/N) for j in N:-1:0];
p = fit(Polynomial, xs, f.(xs));
q = fit(ArnoldiFit, xs, f.(xs));
maximum(abs, p(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 3.304586010148457e16
maximum(abs, q(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 1.1939520722092922e-7
```

```
N = 250; xs = [cos(j*pi/N) for j in N:-1:0];
p = fit(Polynomial, xs, f.(xs));
q = fit(ArnoldiFit, xs, f.(xs));
maximum(abs, p(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 3.55318186254542e92
maximum(abs, q(x) - f(x) for x ∈ range(-1,stop=1,length=500)) # 8.881784197001252e-16
```

```
p = fit(Polynomial, xs, f.(xs), 10); # least-squares fit
q = fit(ArnoldiFit, xs, f.(xs), 10);
maximum(abs, q(x) - p(x) for x ∈ range(-1,stop=1,length=500)) # 4.6775083806238626e-14
Polynomials.norm(q-p, Inf) # 2.2168933355715126e-12 # promotes `q` to `Polynomial`
```

To manipulate the fitted polynomial, conversion is necessary. Conversion can lead to wildly divergent polynomials when n is large.

"""
function polyfitA(x, y, n=length(x)-1; var=:x)
    m = length(x)
    T = eltype(y)
    Q = ones(T, m, n+1)
    #H = UpperHessenberg(zeros(T, n+1, n))
    H = zeros(T, n+1, n)

    q = zeros(T, m)

    # we have Vₓ = QR, y = Vₓa = Q(Ra) = Qd, so d = Q \ y
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

# from Vₓ = QR, we get Vₛ = WR and f = Vₛa = WRa = W(d) stored above
# this finds W
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

This is useful for polynomial evaluation, but other polynomial operations are not defined. Though these fitted polynomials may be converted to other types, for larger degrees this will prove unstable.

"""
struct ArnoldiFit{T, M<:AbstractArray{T,2}, X}  <: AbstractPolynomial{T,X}
    coeffs::Vector{T}
    H::M
end
export ArnoldiFit
@register ArnoldiFit
domain(::Type{<:ArnoldiFit}) = Interval{Open,Open}(-Inf, Inf)

Base.show(io::IO, mimetype::MIME"text/plain", p::ArnoldiFit) = print(io, "ArnoldiFit of degree $(length(p.coeffs)-1)")

evalpoly(x, p::ArnoldiFit) = polyvalA(p.coeffs, p.H, x)

fit(::Type{ArnoldiFit}, x::AbstractVector{T}, y::AbstractVector{T}, deg::Int=length(x)-1;  var=:x, kwargs...) where{T} = polyfitA(x, y, deg; var=var)

Base.convert(::Type{P}, p::ArnoldiFit{T,M,X}) where {P <: AbstractPolynomial,T,M,X} = p(variable(P,indeterminate(p)))




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
# |p(x) - p̃(x)|/|p(x)| ≤ α(n)⋅u ⋅ cond(p,x), where u = finite precision of computation (2^-p)
function LinearAlgebra.cond(p::P, x) where {P <: StandardBasisPolynomial}
    p̃ = map(abs, p)
    p̃(abs(x))/ abs(p(x))
end
