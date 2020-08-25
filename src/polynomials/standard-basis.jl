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
    ⟒(P)(as, p.var)
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
    return ⟒(P)(a2, p.var)
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
    return ⟒(P)(a2, p.var)
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

    return ⟒(P)(q_coeff, var), ⟒(P)(r_coeff, var)

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
function  gcd_noda_sasaki(p::P, q::Q;
                          atol::Real=zero(real(promote_type(T,S))),
                          rtol::Real=Base.rtoldefault(real(promote_type(T,S)))
                          ) where {T,S,
                                   P<: StandardBasisPolynomial{T},
                                   Q<: StandardBasisPolynomial{S},
                                   }
    ⟒(P) == ⟒(Q) ||  return gcd_noda_sasaki(promote(p,q);  atol=atol, rtol=rtol)
    ## check symbol
    a, b = coeffs(p), coeffs(q)
    as =  _gcd_noda_sasaki(a,b, atol=atol,  rtol=rtol)

    ⟒(P)(as, p.var)
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
        return []
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
