export LaurentPolynomial

"""
    LaurentPolynomial(coeffs, range, var)

A [Laurent](https://en.wikipedia.org/wiki/Laurent_polynomial) polynomial is of the form `a_{m}x^m + ... + a_{n}x^n` where `m,n` are  integers (not necessarily positive) with ` m <= n`.

The `coeffs` specify `a_{m}, a_{m-1}, ..., a_{n}`. Rhe range specified is of the  form  `m:n`,  if left  empty, `0:length(coeffs)-1` is  used (i.e.,  the coefficients refer  to the standard basis).

Laurent polynomials and standard basis polynomials  promote to  Laurent polynomials. Laurent polynomials may be  converted to a standard basis  polynomial when `m >= 0`
. 

Integration will fail if there is a `x⁻¹` term in the polynomial.

Example:
```jldoctest
julia> using Polynomials

julia> P = LaurentPolynomial
LaurentPolynomial

julia> p = P([1,1,1],  -1:1)
LaurentPolynomial(x⁻¹ + 1 + x)

julia> q = P([1,1,1])
LaurentPolynomial(1 + x + x²)

julia> pp = Polynomial([1,1,1])
Polynomial(1 + x + x^2)

julia> p + q
LaurentPolynomial(x⁻¹ + 2 + 2*x + x²)

julia> p * q
LaurentPolynomial(x⁻¹ + 2 + 3*x + 2*x² + x³)

julia> p * pp
LaurentPolynomial(x⁻¹ + 2 + 3*x + 2*x² + x³)

julia> pp - q
LaurentPolynomial(0)

julia> derivative(p)
LaurentPolynomial(-x⁻² + 1)

julia> integrate(q)
LaurentPolynomial(1.0*x + 0.5*x² + 0.3333333333333333*x³)

julia> integrate(p)  # x⁻¹  term is an issue
ERROR: ArgumentError: Can't integrate Laurent  polynomial with  `x⁻¹` term

julia> integrate(P([1,1,1], -5:-3))
LaurentPolynomial(-0.25*x⁻⁴ - 0.3333333333333333*x⁻³ - 0.5*x⁻²)

julia> x⁻¹ = inv(variable(LaurentPolynomial)) # `inv` defined on monomials
LaurentPolynomial(1.0*x⁻¹)

julia> p = Polynomial([1,2,3])
Polynomial(1 + 2*x + 3*x^2)

julia> x = variable()
Polynomial(x)

julia> x^degree(p) * p(x⁻¹) # reverses  coefficients
LaurentPolynomial(3.0 + 2.0*x + 1.0*x²)
```
"""
struct LaurentPolynomial{T <: Number} <: StandardBasisPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    m::Base.RefValue{Int64}
    n::Base.RefValue{Int64}
    function LaurentPolynomial{T}(coeffs::AbstractVector{T},
                                  rng::UnitRange{Int64}=0:length(coeffs)-1,
                                  var::Symbol=:x) where {T <: Number}

        m,n = first(rng), last(rng)

        # trim zeros from front and back
        lnz = findlast(!iszero, coeffs)
        fnz = findfirst(!iszero, coeffs)
        (lnz == nothing || length(coeffs) == 0) && return new{T}(zeros(T,1), var, Ref(0), Ref(0))
        if  lnz !=  length(rng) ||  fnz !=  1
            coeffs =  coeffs[fnz:lnz]
            m = m + fnz - 1
            n = m + (lnz-fnz)
        end
        (n-m+1  == length(coeffs)) || throw(ArgumentError("Lengths do not match"))

        new{T}(coeffs, var, Ref(m),  Ref(n))

    end
end

@register LaurentPolynomial

function  LaurentPolynomial{T}(coeffs::AbstractVector{S},
                               rng::UnitRange{Int64}=0:length(coeffs)-1,
                               var::Symbol=:x) where {T <: Number, S <: Number}
    LaurentPolynomial{T}(T.(coeffs), rng, var)
end

function LaurentPolynomial(coeffs::AbstractVector{T}, rng::UnitRange, var::SymbolLike=:x) where {T <: Number}
    LaurentPolynomial{T}(coeffs, rng, Symbol(var))
end

function LaurentPolynomial(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {T <: Number}
    LaurentPolynomial{T}(coeffs, 0:length(coeffs)-1, Symbol(var))
end

## Alternate interface
Polynomial(coeffs::AbstractVector{T}, rng::UnitRange, var::SymbolLike=:x) where {T <: Number} =
    LaurentPolynomial{T}(coeffs, rng, Symbol(var))

Polynomial{T}(coeffs::AbstractVector{S}, rng::UnitRange, var::SymbolLike=:x) where {T <: Number, S <: Number} =
    LaurentPolynomial{T}(T.(coeffs), rng, Symbol(var))

##
## conversion
##

# LaurentPolynomial is a wider collection than other standard basis polynomials.
Base.promote_rule(::Type{P},::Type{Q}) where {T, P <: LaurentPolynomial{T}, S, Q <: StandardBasisPolynomial{S}} =
    LaurentPolynomial{promote_type(T, S)}

Base.promote_rule(::Type{Q},::Type{P}) where {T, P <: LaurentPolynomial{T}, S, Q <: StandardBasisPolynomial{S}} =
    LaurentPolynomial{promote_type(T, S)}

function Base.convert(P::Type{<:Polynomial}, q::LaurentPolynomial)
    m,n = extrema(q)
    m < 0 && throw(ArgumentError("Can't convert a Laurent polynomial with m < 0"))
    P([q[i] for i  in 0:n], q.var)
end

function Base.convert(::Type{P}, q::StandardBasisPolynomial{S}) where {T, P <:LaurentPolynomial{T},S}
    d = degree(q)
    P([q[i] for i in 0:d], 0:degree(q), q.var)
end

##
## generic functions
##
Base.extrema(p::LaurentPolynomial) =  (p.m[],p.n[])
Base.range(p::LaurentPolynomial) = p.m[]:p.n[]
function Base.inv(p::LaurentPolynomial)
    m,n =  extrema(p)
    m != n && throw(ArgumentError("Only monomials can be inverted"))
    LaurentPolynomial([1/p for p in p.coeffs], -m:-n, p.var)
end

##
## changes to common.jl mostly as the range in the type is different
##
Base.copy(p::P) where {P <: LaurentPolynomial} = P(copy(coeffs(p)), range(p), p.var)
Base.:(==)(p1::LaurentPolynomial, p2::LaurentPolynomial) =
    check_same_variable(p1, p2) && (range(p1) == range(p2)) && (coeffs(p1) == coeffs(p2))
Base.hash(p::LaurentPolynomial, h::UInt) = hash(p.var, hash(range(p), hash(coeffs(p), h)))

degree(p::LaurentPolynomial) = p.n[]
isconstant(p::LaurentPolynomial) = range(p) == 0:0
basis(P::Type{<:LaurentPolynomial{T}}, n::Int, var::SymbolLike=:x) where{T} = LaurentPolynomial(ones(T,1), n:n, var)
basis(P::Type{LaurentPolynomial}, n::Int, var::SymbolLike=:x) = LaurentPolynomial(ones(Float64, 1), n:n, var)

Base.zero(::Type{LaurentPolynomial{T}},  var=Symbollike=:x) where {T} =  LaurentPolynomial{T}(zeros(T,1),  0:0, Symbol(var))
Base.zero(::Type{LaurentPolynomial},  var=Symbollike=:x) =  zero(LaurentPolynomial{Float64}, var)
Base.zero(p::P, var=Symbollike=:x) where {P  <: LaurentPolynomial} = zero(P, var)


# get/set index. Work with  offset
function Base.getindex(p::LaurentPolynomial{T}, idx::Int) where {T <: Number}
    m,n = extrema(p)
    i = idx - m + 1
    (i < 1 || i > (n-m+1))  && return zero(T)
    p.coeffs[i]
end

# extend if out of bounds
function Base.setindex!(p::LaurentPolynomial{T}, value::Number, idx::Int) where {T}

    m,n = extrema(p)
    if idx > n
        append!(p.coeffs, zeros(T, idx-n))
        n = idx
        p.n[] = n
    elseif idx < m
        prepend!(p.coeffs, zeros(T, m-idx))
        m = idx
        p.m[] = m
    end

    i = idx - m + 1
    p.coeffs[i] = value

    return p
    
end

Base.firstindex(p::LaurentPolynomial) = p.m[]
Base.lastindex(p::LaurentPolynomial) = p.n[]
Base.eachindex(p::LaurentPolynomial) = range(p)

## chop/truncation
# trim  from *both* ends
function chop!(p::P;
               rtol::Real = Base.rtoldefault(real(T)),
               atol::Real = 0,) where {T, P <: LaurentPolynomial{T}}
    
    m0,n0 = m,n = extrema(p)
    for k in n:-1:m
        if isapprox(p[k], zero(T); rtol = rtol, atol = atol)
            n -= 1
        else
            break
        end
    end
    for k in m:n-1
        if isapprox(p[k], zero(T); rtol = rtol, atol = atol)
            m += 1
        else
            break
        end
    end

    cs = coeffs(p)
    rng = m-m0+1:n-m0+1
    resize!(p.coeffs, length(rng))
    p.coeffs[:] = coeffs(p)[rng]
    isempty(p.coeffs) && push!(p.coeffs,zero(T))
    p.m[], p.n[] = m, max(m,n)

    p

end

function truncate!(p::LaurentPolynomial{T};
                  rtol::Real = Base.rtoldefault(real(T)),
                  atol::Real = 0,) where {T}
    
    max_coeff = maximum(abs, coeffs(p))
    thresh = max_coeff * rtol + atol

    for i in eachindex(p)
        if abs(p[i]) <= thresh
            p[i] = zero(T)
        end
    end

    chop!(p)

    return p
    
end

# use unicode exponents. XXX modify printexponent to always use these?
function showterm(io::IO, ::Type{<:LaurentPolynomial}, pj::T, var, j, first::Bool, mimetype) where {T}
    if iszero(pj) return false end
    pj = printsign(io, pj, first, mimetype)
    if !(pj == one(T) && !(showone(T) || j == 0))
        printcoefficient(io, pj, j, mimetype)
    end
    printproductsign(io, pj, j, mimetype)
    iszero(j) && return  true
    print(io, var)
    j ==1  && return  true
    unicode_exponent(io, j)
    return true
end


##
## ----
##


# evaluation uses `evalpoly`
function (p::LaurentPolynomial{T})(x::S) where {T,S}
    m,n = extrema(p)
    m  == n == 0 && return p[0] * _one(S)
    if m >= 0
        evalpoly(x, NTuple{n+1,T}(p[i] for i in 0:n))
    elseif n <= 0
        evalpoly(inv(x), NTuple{m+1,T}(p[i] for i in 0:-1:m))
    else
        # eval pl(x) = a_mx^m + ...+ a_0 at 1/x; pr(x) = a_0 + a_1x + ... + a_nx^n  at  x; subtract a_0
        evalpoly(inv(x), NTuple{-m+1,T}(p[i] for i in 0:-1:m)) + evalpoly(x, NTuple{n+1,T}(p[i] for i in 0:n)) - p[0]
    end
end
                 
        

# scalar operattoinis
Base.:-(p::P) where {P <: LaurentPolynomial} = P(-coeffs(p), range(p), p.var)

function Base.:+(p::LaurentPolynomial{T}, c::S) where {T, S <: Number}
    q = LaurentPolynomial([c], 0:0, p.var)
    p + q
end

function Base.:*(p::P, c::S) where {T,P <: LaurentPolynomial,  S <: Number}
    as = c * copy(coeffs(p))
    return ⟒(P)(as, range(p), p.var)
end


function Base.:/(p::P, c::S) where {T,P <: LaurentPolynomial{T},S <: Number}
    R = promote_type(P, eltype(one(T) / one(S)))
    return R(coeffs(p) ./ c, range(p), p.var)
end

##
## Poly + and  *
##
function Base.:+(p1::P1, p2::P2) where {T,P1<:LaurentPolynomial{T}, S, P2<:LaurentPolynomial{S}}

    if isconstant(p1)
        p1 = P1(p1.coeffs, range(p1), p2.var)
    elseif isconstant(p2)
        p2 = P2(p2.coeffs, range(p2), p1.var)
    end
    
    p1.var != p2.var && error("LaurentPolynomials must have same variable")

    R = promote_type(T,S)

    m1,n1 = extrema(p1)
    m2,n2 = extrema(p2)
    m,n = min(m1,m2), max(n1, n2)

    as = zeros(R, length(m:n))
    for i in m:n
        as[1 + i-m] = p1[i] + p2[i]
    end

    q = LaurentPolynomial{R}(as, m:n, p1.var)
    chop!(q)

    return q
    
end

function Base.:*(p1::LaurentPolynomial{T}, p2::LaurentPolynomial{S}) where {T,S}

    isconstant(p1) && return p2 * p1[0]
    isconstant(p2) && return p1 * p2[0]

    p1.var != p2.var && error("LaurentPolynomials must have same variable")

    R = promote_type(T,S)

    m1,n1 = extrema(p1)
    m2,n2 = extrema(p2)
    m,n = m1 + m2, n1+n2

    as = zeros(R, length(m:n))
    for i in eachindex(p1)
        p1ᵢ = p1[i]
        for j in eachindex(p2)
            as[1 + i+j - m] = muladd(p1ᵢ, p2[j], as[1 + i + j - m])
        end
    end

    p = LaurentPolynomial(as, m:n, p1.var)
    chop!(p)

    return p
end

##
## d/dx, ∫
##
function derivative(p::P, order::Integer = 1) where {T, P<:LaurentPolynomial{T}}
    
    order < 0 && error("Order of derivative must be non-negative")
    order == 0 && return p

    hasnan(p) && return ⟒(P)(T[NaN], 0:0, p.var)

    m,n = extrema(p)
    m = m - order
    n = n - order
    as =  zeros(T, length(m:n))

    for k in eachindex(p)
        idx = 1 + k - order - m
        if 0 ≤ k ≤ order - 1
            as[idx] = zero(T)
        else
            as[idx] = reduce(*, (k - order + 1):k, init = p[k])
        end
    end
    
    chop!(LaurentPolynomial(as, m:n, p.var))
    
end


function integrate(p::P, k::S) where {T, P<: LaurentPolynomial{T}, S<:Number}

    !iszero(p[-1])  && throw(ArgumentError("Can't integrate Laurent  polynomial with  `x⁻¹` term"))
    R = eltype((one(T)+one(S))/1)

    if hasnan(p) || isnan(k)
        return P([NaN], 0:0, p.var) # not R(NaN)!! don't like XXX
    end


    m,n = extrema(p)
    if  n < 0
        n = 0
    else
        n += 1
    end
    if m < 0
        m += 1
    else
        m = 0
    end
    as = zeros(R,  length(m:n))
    
    for k in eachindex(p)
        as[1 + k+1-m]  =  p[k]/(k+1)
    end

    as[1-m] = k

    return ⟒(P)(as, m:n, p.var)
    
end


function Base.gcd(p::LaurentPolynomial{T}, q::LaurentPolynomial{T}, args...; kwargs...) where {T}
    mp, Mp = extrema(p)
    mq, Mq = extrema(q)
    if mp < 0 || mq < 0
        throw(ArgumentError("GCD is not defined when there are `x⁻¹` terms"))
    end

    degree(p) == 0 && return iszero(p) ? q : one(q)
    degree(q) == 0 && return iszero(q) ? p : one(p)
    check_same_variable(p,q) || throw(ArgumentError("p and q have different symbols"))
    
    pp, qq = convert(Polynomial, p), convert(Polynomial, q)
    u = gcd(pp, qq, args..., kwargs...)
    return LaurentPolynomial(coeffs(u), p.var)
end
