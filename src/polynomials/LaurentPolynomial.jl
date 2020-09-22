export LaurentPolynomial

"""
    LaurentPolynomial(coeffs, range, var)

A [Laurent](https://en.wikipedia.org/wiki/Laurent_polynomial) polynomial is of the form `a_{m}x^m + ... + a_{n}x^n` where `m,n` are  integers (not necessarily positive) with ` m <= n`.

The `coeffs` specify `a_{m}, a_{m-1}, ..., a_{n}`. The range specified is of the  form  `m`,  if left  empty, `m` is taken to be `0` (i.e.,  the coefficients refer  to the standard basis). Alternatively, the coefficients can be specified using an `OffsetVector` from the `OffsetArrays` package.

Laurent polynomials and standard basis polynomials  promote to  Laurent polynomials. Laurent polynomials may be  converted to a standard basis  polynomial when `m >= 0`
.

Integration will fail if there is a `x⁻¹` term in the polynomial.

Example:
```jldoctest
julia> using Polynomials

julia> P = LaurentPolynomial
LaurentPolynomial

julia> p = P([1,1,1],  -1)
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

julia> integrate(P([1,1,1], -5))
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
                                  m::Int,
                                  var::Symbol=:x) where {T <: Number}

        
        # trim zeros from front and back
        lnz = findlast(!iszero, coeffs)
        fnz = findfirst(!iszero, coeffs)
        (lnz == nothing || length(coeffs) == 0) && return new{T}(zeros(T,1), var, Ref(0), Ref(0))
        coeffs =  coeffs[fnz:lnz]
        m = m + fnz - 1
        n = m + (lnz-fnz)

        (n-m+1  == length(coeffs)) || throw(ArgumentError("Lengths do not match"))

        new{T}(coeffs, var, Ref(m),  Ref(n))

    end
end

@register LaurentPolynomial

## constructors
function LaurentPolynomial{T}(coeffs::AbstractVector{S}, m::Int, var::SymbolLike=:x) where {
    T <: Number, S <: Number}
    LaurentPolynomial{T}(T.(coeffs), m, var)
end

function LaurentPolynomial{T}(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {
    T <: Number}
    LaurentPolynomial{T}(coeffs, 0, var)
end

function  LaurentPolynomial{T}(coeffs::OffsetArray{T, 1, Array{T,1}}, var::SymbolLike=:x) where {
    T<:Number}
    m,n = axes(coeffs, 1)
    LaurentPolynomial{T}(T.(coeffs.parent), m, Symbol(var))
end

function LaurentPolynomial(coeffs::AbstractVector{T}, m::Int, var::SymbolLike=:x) where {T <: Number}
    LaurentPolynomial{T}(coeffs, m, Symbol(var))
end




## Alternate with range specified
## Deprecate
function  LaurentPolynomial{T}(coeffs::AbstractVector{S},
                               rng::UnitRange{Int64},
                               var::Symbol=:x) where {T <: Number, S <: Number}
    Base.depwarn("Using a range to indicate the offset is deprecated. Use just the lower value",
                 :LaurentPolynomial)
    error("")
    LaurentPolynomial{T}(T.(coeffs), first(rng), var)
end

function LaurentPolynomial(coeffs::AbstractVector{T}, rng::UnitRange, var::SymbolLike=:x) where {T <: Number}
    Base.depwarn("Using a range to indicate the offset is deprecated. Use just the lower value",
                 :LaurentPolynomial)
    LaurentPolynomial{T}(coeffs, rng, Symbol(var))
end



##
## conversion
##

# LaurentPolynomial is a wider collection than other standard basis polynomials.
Base.promote_rule(::Type{P},::Type{Q}) where {T, P <: LaurentPolynomial{T}, S, Q <: StandardBasisPolynomial{S}} =
    LaurentPolynomial{promote_type(T, S)}

Base.promote_rule(::Type{Q},::Type{P}) where {T, P <: LaurentPolynomial{T}, S, Q <: StandardBasisPolynomial{S}} =
    LaurentPolynomial{promote_type(T, S)}

function Base.convert(P::Type{<:Polynomial}, q::LaurentPolynomial)
    m,n = (extrema∘degreerange)(q)
    m < 0 && throw(ArgumentError("Can't convert a Laurent polynomial with m < 0"))
    P([q[i] for i  in 0:n], q.var)
end

function Base.convert(::Type{P}, q::StandardBasisPolynomial{S}) where {T, P <:LaurentPolynomial{T},S}
    d = degree(q)
    P([q[i] for i in 0:d], 0, q.var)
end

##
## generic functions
##
function Base.extrema(p::LaurentPolynomial)
    Base.depwarn("`extrema(::LaurentPolynomial)` is deprecated. Use `(firstindex(p), lastindex(p))`", :extrema)
    (p.m[], p.n[])
end
function Base.range(p::LaurentPolynomial)
    Base.depwarn("`range(::LaurentPolynomial)` is deprecated. Use `firstindex(p):lastindex(p)`", :range)
    p.m[]:p.n[]
end

function Base.inv(p::LaurentPolynomial)
    m,n =  (extrema∘degreerange)(p)
    m != n && throw(ArgumentError("Only monomials can be inverted"))
    LaurentPolynomial([1/p for p in p.coeffs], -m, p.var)
end

##
## changes to common.jl mostly as the range in the type is different
##
Base.:(==)(p1::LaurentPolynomial, p2::LaurentPolynomial) =
    check_same_variable(p1, p2) && (degreerange(p1) == degreerange(p2)) && (coeffs(p1) == coeffs(p2))
Base.hash(p::LaurentPolynomial, h::UInt) = hash(p.var, hash(degreerange(p), hash(coeffs(p), h)))

isconstant(p::LaurentPolynomial) = iszero(lastindex(p)) && iszero(firstindex(p)) 
basis(P::Type{<:LaurentPolynomial{T}}, n::Int, var::SymbolLike=:x) where{T} = LaurentPolynomial(ones(T,1), n, var)
basis(P::Type{LaurentPolynomial}, n::Int, var::SymbolLike=:x) = LaurentPolynomial(ones(Float64, 1), n, var)

Base.zero(::Type{LaurentPolynomial{T}},  var=Symbollike=:x) where {T} =  LaurentPolynomial{T}(zeros(T,1),  0, Symbol(var))
Base.zero(::Type{LaurentPolynomial},  var=Symbollike=:x) =  zero(LaurentPolynomial{Float64}, var)
Base.zero(p::P, var=Symbollike=:x) where {P  <: LaurentPolynomial} = zero(P, var)


# get/set index. Work with  offset
function Base.getindex(p::LaurentPolynomial{T}, idx::Int) where {T <: Number}
    m,n = (extrema ∘ degreerange)(p)
    i = idx - m + 1
    (i < 1 || i > (n-m+1))  && return zero(T)
    p.coeffs[i]
end

# extend if out of bounds
function Base.setindex!(p::LaurentPolynomial{T}, value::Number, idx::Int) where {T}

    m,n = (extrema ∘ degreerange)(p)
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
Base.eachindex(p::LaurentPolynomial) = degreerange(p)
degreerange(p::LaurentPolynomial) = firstindex(p):lastindex(p)

_convert(p::P, as) where {P <: LaurentPolynomial} = ⟒(P)(as, firstindex(p), p.var)

## chop/truncation
# trim  from *both* ends
function chop!(p::P;
               rtol::Real = Base.rtoldefault(real(T)),
               atol::Real = 0,) where {T, P <: LaurentPolynomial{T}}

    m0,n0 = m,n = (extrema ∘ degreerange)(p)
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
## ---- Conjugation has different defintions
##

"""
    conj(p)

This satisfies `conj(p(x)) = conj(p)(conj(x)) = p̄(conj(x))` or `p̄(x) = (conj ∘ p ∘ conj)(x)`

Examples
```jldoctest
julia> z = variable(LaurentPolynomial, :z)
LaurentPolynomial(z)

julia> p = LaurentPolynomial([im, 1+im, 2 + im], -1:1, :z)
LaurentPolynomial(im*z⁻¹ + (1 + 1im) + (2 + 1im)*z)

julia> conj(p)(conj(z)) ≈ conj(p(z))
true

julia> conj(p)(z) ≈ (conj ∘ p ∘ conj)(z)
true
```
"""
LinearAlgebra.conj(p::P) where {P <: LaurentPolynomial} = map(conj, p)


"""
    paraconj(p)

[cf.](https://ccrma.stanford.edu/~jos/filters/Paraunitary_FiltersC_3.html)

Call `p̂ = paraconj(p)` and `p̄` = conj(p)`, then this satisfies
`conj(p(z)) = p̂(1/conj(z))` or `p̂(z) = p̄(1/z) = (conj ∘ p ∘ conj ∘ inf)(z)`.

Examples:

```jldoctest
julia> z = variable(LaurentPolynomial, :z)
LaurentPolynomial(z)

julia> h = LaurentPolynomial([1,1], -1:0, :z)
LaurentPolynomial(z⁻¹ + 1)

julia> Polynomials.paraconj(h)(z) ≈ 1 + z ≈ LaurentPolynomial([1,1], 0:1, :z)
true

julia> h = LaurentPolynomial([3,2im,1], -2:0, :z)
LaurentPolynomial(3*z⁻² + 2im*z⁻¹ + 1)

julia> Polynomials.paraconj(h)(z) ≈ 1 - 2im*z + 3z^2 ≈ LaurentPolynomial([1, -2im, 3], 0:2, :z)
true

julia> Polynomials.paraconj(h)(z) ≈ (conj ∘ h ∘ conj ∘ inv)(z)
true
"""
function paraconj(p::LaurentPolynomial)
    cs = p.coeffs
    ds = adjoint.(cs)
    n = degree(p)
    LaurentPolynomial(reverse(ds), -n, p.var)
end

"""
    cconj(p)

Conjugation of a polynomial with respect to the imaginary axis.

The `cconj` of a polynomial, `p̃`, conjugates the coefficients and applies `s -> -s`. That is `cconj(p)(s) = conj(p)(-s)`.

This satisfies for *imaginary* `s`: `conj(p(s)) = p̃(s) = (conj ∘ p)(s) = cconj(p)(s) `

[ref](https://github.com/hurak/PolynomialEquations.jl#symmetrix-conjugate-equation-continuous-time-case)

Examples:
```jldoctest
julia> s = 2im
0 + 2im

julia> p = LaurentPolynomial([im,-1, -im, 1], 1:2, :s)
LaurentPolynomial(im*s - s² - im*s³ + s⁴)

julia> Polynomials.cconj(p)(s) ≈ conj(p(s))
true

julia> a = LaurentPolynomial([-0.12, -0.29, 1],:s)
LaurentPolynomial(-0.12 - 0.29*s + 1.0*s²)

julia> b = LaurentPolynomial([1.86, -0.34, -1.14, -0.21, 1.19, -1.12],:s)
LaurentPolynomial(1.86 - 0.34*s - 1.14*s² - 0.21*s³ + 1.19*s⁴ - 1.12*s⁵)

julia> x = LaurentPolynomial([-15.5, 50.0096551724139, 1.19], :s)
LaurentPolynomial(-15.5 + 50.0096551724139*s + 1.19*s²)

julia> Polynomials.cconj(a) * x + a * Polynomials.cconj(x) ≈ b + Polynomials.cconj(b)
true
```

"""
function cconj(p::LaurentPolynomial)
    ps = conj.(coeffs(p))
    m,n = (extrema ∘ degreerange)(p)
    for i in m:n
        if isodd(i)
            ps[i+1-m] *= -1
        end
    end
    LaurentPolynomial(ps, m, p.var)
end



##
## ----
##


# evaluation uses `evalpoly`
function (p::LaurentPolynomial{T})(x::S) where {T,S}
    m,n = (extrema ∘ degreerange)(p)
    m  == n == 0 && return p[0] * _one(S)
    if m >= 0
        evalpoly(x, NTuple{n+1,T}(p[i] for i in 0:n))
    elseif n <= 0
        evalpoly(inv(x), NTuple{-m+1,T}(p[i] for i in 0:-1:m))
    else
        # eval pl(x) = a_mx^m + ...+ a_0 at 1/x; pr(x) = a_0 + a_1x + ... + a_nx^n  at  x; subtract a_0
        l = evalpoly(inv(x), NTuple{-m+1,T}(p[i] for i in 0:-1:m))
        r =  evalpoly(x, NTuple{n+1,T}(p[i] for i in 0:n))
        mid = p[0]
        l + r - mid
    end
end



# scalar operattoinis
# standard-basis defn. assumes basis 1, x, x², ...
Base.:+(p::LaurentPolynomial{T}, c::S) where {T, S <: Number} = sum(promote(p,c))

##
## Poly + and  *
##
function Base.:+(p1::P1, p2::P2) where {T,P1<:LaurentPolynomial{T}, S, P2<:LaurentPolynomial{S}}

    if isconstant(p1)
        p1 = P1(p1.coeffs, firstindex(p1), p2.var)
    elseif isconstant(p2)
        p2 = P2(p2.coeffs, firstindex(p2), p1.var)
    end

    p1.var != p2.var && error("LaurentPolynomials must have same variable")

    R = promote_type(T,S)

    m1,n1 = (extrema ∘ degreerange)(p1)
    m2,n2 = (extrema ∘ degreerange)(p2)
    m,n = min(m1,m2), max(n1, n2)

    as = zeros(R, length(m:n))
    for i in m:n
        as[1 + i-m] = p1[i] + p2[i]
    end

    q = LaurentPolynomial{R}(as, m, p1.var)
    chop!(q)

    return q

end

function Base.:*(p1::LaurentPolynomial{T}, p2::LaurentPolynomial{S}) where {T,S}

    isconstant(p1) && return p2 * p1[0]
    isconstant(p2) && return p1 * p2[0]

    p1.var != p2.var && error("LaurentPolynomials must have same variable")

    R = promote_type(T,S)

    m1,n1 = (extrema ∘ degreerange)(p1)
    m2,n2 = (extrema ∘ degreerange)(p2)
    m,n = m1 + m2, n1+n2

    as = zeros(R, length(m:n))
    for i in eachindex(p1)
        p1ᵢ = p1[i]
        for j in eachindex(p2)
            as[1 + i+j - m] = muladd(p1ᵢ, p2[j], as[1 + i + j - m])
        end
    end

    p = LaurentPolynomial(as, m, p1.var)
    chop!(p)

    return p
end

##
## roots
##
"""
    roots(p)

Compute the roots of the Laurent polynomial `p`.

The roots of a function (Laurent polynomial in this case) `a(z)` are the values of `z` for which the function vanishes. A Laurent polynomial ``a(z) = a_m z^m + a_{m+1} z^{m+1} + ... + a_{-1} z^{-1} + a_0 + a_1 z + ... + a_{n-1} z^{n-1} + a_n z^n`` can equivalently be viewed as a rational function with a multiple singularity (pole) at the origin. The roots are then the roots of the numerator polynomial. For example, ``a(z) = 1/z + 2 + z`` can be written as ``a(z) = (1+2z+z^2) / z`` and the roots of `a` are the roots of ``1+2z+z^2``.

# Example

```julia
julia> p = LaurentPolynomial([24,10,-15,0,1],-2:1,:z)
LaurentPolynomial(24*z⁻² + 10*z⁻¹ - 15 + z²)

julia> roots(a)
4-element Array{Float64,1}:
 -3.999999999999999
 -0.9999999999999994
  1.9999999999999998
  2.9999999999999982
```
"""
function  roots(p::P; kwargs...)  where  {T, P <: LaurentPolynomial{T}}
    c = coeffs(p)
    r = degreerange(p)
    d = r[end] - min(0, r[1]) + 1    # Length of the coefficient vector, taking into consideration
                                     # the case when the lower degree is strictly positive
                                     # (like p=3z^2).
    z = zeros(T, d)                  # Reserves space for the coefficient vector.
    z[max(0, r[1]) + 1:end] = c      # Leaves the coeffs of the lower powers as zeros.
    a = Polynomial(z, p.var)         # The root is then the root of the numerator polynomial.
    return roots(a; kwargs...)
end

##
## d/dx, ∫
##
function derivative(p::P, order::Integer = 1) where {T, P<:LaurentPolynomial{T}}

    order < 0 && error("Order of derivative must be non-negative")
    order == 0 && return p

    hasnan(p) && return ⟒(P)(T[NaN], 0, p.var)

    m,n = (extrema ∘ degreerange)(p)
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

    chop!(LaurentPolynomial(as, m, p.var))

end


function integrate(p::P, k::S) where {T, P<: LaurentPolynomial{T}, S<:Number}

    !iszero(p[-1])  && throw(ArgumentError("Can't integrate Laurent  polynomial with  `x⁻¹` term"))
    R = eltype((one(T)+one(S))/1)

    if hasnan(p) || isnan(k)
        return P([NaN], 0, p.var) # not R(NaN)!! don't like XXX
    end


    m,n = (extrema ∘ degreerange)(p)
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

    return ⟒(P)(as, m, p.var)

end


function Base.gcd(p::LaurentPolynomial{T}, q::LaurentPolynomial{T}, args...; kwargs...) where {T}
    mp, Mp = (extrema ∘ degreerange)(p)
    mq, Mq = (extrema ∘ degreerange)(q)
    if mp < 0 || mq < 0
        throw(ArgumentError("GCD is not defined when there are `x⁻ⁿ` terms"))
    end

    degree(p) == 0 && return iszero(p) ? q : one(q)
    degree(q) == 0 && return iszero(q) ? p : one(p)
    check_same_variable(p,q) || throw(ArgumentError("p and q have different symbols"))

    pp, qq = convert(Polynomial, p), convert(Polynomial, q)
    u = gcd(pp, qq, args..., kwargs...)
    return LaurentPolynomial(coeffs(u), p.var)
end
