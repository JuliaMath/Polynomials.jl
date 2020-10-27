abstract type StandardBasisPolynomial{T} <: AbstractPolynomial{T} end



function showterm(io::IO, ::Type{<:StandardBasisPolynomial}, pj::T, var, j, first::Bool, mimetype) where {T}

    if iszero(pj) return false end

    pj = printsign(io, pj, first, mimetype)

    if !(pj == one(T) && !(showone(T) || j == 0))
        printcoefficient(io, pj, j, mimetype)
    end

    printproductsign(io, pj, j, mimetype)
    printexponent(io, var, j, mimetype)
    return true
end

# allows  broadcast  issue #209
evalpoly(x, p::StandardBasisPolynomial) = p(x)

domain(::Type{<:StandardBasisPolynomial}) = Interval(-Inf, Inf)
mapdomain(::Type{<:StandardBasisPolynomial}, x::AbstractArray) = x

## generic test if polynomial `p` is a constant
isconstant(p::StandardBasisPolynomial) = degree(p) <= 0

Base.convert(P::Type{<:StandardBasisPolynomial}, q::StandardBasisPolynomial) = isa(q, P) ? q : P([q[i] for i in 0:degree(q)], q.var)

Base.values(p::StandardBasisPolynomial) = values(p.coeffs)

variable(::Type{P}, var::SymbolLike = :x) where {P <: StandardBasisPolynomial} = P([0, 1], var)

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


function Base.:+(p::P, c::S) where {T, P <: StandardBasisPolynomial{T}, S<:Number}
    R = promote_type(T,S)
    as = R[c  for c in coeffs(p)]
    as[1] += c
    _convert(p, as)
end


function derivative(p::P, order::Integer = 1) where {T, P <: StandardBasisPolynomial{T}}
    order < 0 && error("Order of derivative must be non-negative")

    # we avoid usage like Base.promote_op(*, T, Int) here, say, as
    # Base.promote_op(*, Rational, Int) is Any, not Rational in analogy to
    # Base.promote_op(*, Complex, Int)
    R = eltype(one(T)*1)
    order == 0 && return p
    hasnan(p) && return ⟒(P){R}(R[NaN], p.var)
    order > length(p) && return zero(⟒(P){R},p.var)
    d = degree(p)
    d <= 0 && return zero(⟒(P){R},p.var)
    n = d + 1
    a2 = Vector{R}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return _convert(p, a2)
end


function integrate(p::P, k::S) where {T, P <: StandardBasisPolynomial{T}, S<:Number}

    R = eltype((one(T)+one(S))/1)
    if hasnan(p) || isnan(k)
        return ⟒(P)([NaN])
    end
    n = length(p)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return _convert(p, a2)
end


function Base.divrem(num::P, den::Q) where {T, P <: StandardBasisPolynomial{T}, S, Q <: StandardBasisPolynomial{S}}

    check_same_variable(num, den) || error("Polynomials must have same variable")
    var = num.var


    n = degree(num)
    m = degree(den)

    m == -1 && throw(DivideError())
    if m == 0 && den[0] ≈ 0 throw(DivideError()) end

    R = eltype(one(T)/one(S))

    deg = n - m + 1

    if deg ≤ 0
        return zero(P, var), num
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
    d < 1 && error("Series must have degree greater than 1")
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
    comp  = companion(⟒(P)(as[k:K], p.var))
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
function compensated_horner(p::P, x) where {T, P <: Polynomials.StandardBasisPolynomial{T}}
    compensated_horner(coeffs(p), x)
end

# rexpressed from paper to compute horner_sum in same pass
# sᵢ -- horner sum
# c -- compensating term
@inline function compensated_horner(ps, x) 
    n, T = length(ps), eltype(ps)
    aᵢ = ps[end]
    sᵢ = aᵢ * _one(x)
    c = zero(T) * _one(x)
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
        sσᵢ =:(ps[end] * _one(x), zero(S))
        c = :(zero(S) * _one(x))
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
function LinearAlgebra.cond(p::P, x) where {P <: Polynomials.StandardBasisPolynomial}
    p̃ = map(abs, p)
    p̃(abs(x))/ abs(p(x))
end



