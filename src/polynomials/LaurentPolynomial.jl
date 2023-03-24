export LaurentPolynomial


"""
    LaurentPolynomial{T,X}(coeffs::AbstractVector, [m::Integer = 0], [var = :x])

A [Laurent](https://en.wikipedia.org/wiki/Laurent_polynomial) polynomial is of the form `a_{m}x^m + ... + a_{n}x^n` where `m,n` are  integers (not necessarily positive) with ` m <= n`.

The `coeffs` specify `a_{m}, a_{m-1}, ..., a_{n}`.
The argument `m` represents the lowest exponent of the variable in the series, and is taken to be zero by default.

Laurent polynomials and standard basis polynomials promote to  Laurent polynomials. Laurent polynomials may be  converted to a standard basis  polynomial when `m >= 0`
.

Integration will fail if there is a `x⁻¹` term in the polynomial.

!!! note
    `LaurentPolynomial` is axis-aware, unlike the other polynomial types in this package.

# Examples:
```jldoctest laurent
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
ERROR: ArgumentError: Can't integrate Laurent polynomial with  `x⁻¹` term

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
struct LaurentPolynomial{T, X} <: LaurentBasisPolynomial{T, X}
    coeffs::Vector{T}
    m::Base.RefValue{Int}
    n::Base.RefValue{Int}
    function LaurentPolynomial{T,X}(coeffs::AbstractVector{S},
                                    m::Union{Int, Nothing}=nothing) where {T, X, S}

        fnz = findfirst(!iszero, coeffs)
        isnothing(fnz) && return new{T,X}(zeros(T,1), Ref(0), Ref(0))
        lnz = findlast(!iszero, coeffs)
        if Base.has_offset_axes(coeffs)
            # if present, use axes
            cs = convert(Vector{T}, coeffs[fnz:lnz])
            return new{T,X}(cs, Ref(fnz), Ref(lnz))
        else

            c = convert(Vector{T}, coeffs[fnz:lnz])

            m′ = fnz - 1 + (isnothing(m) ? 0 : m)
            n = m′ + (lnz-fnz)

            (n - m′ + 1  == length(c)) || throw(ArgumentError("Lengths do not match"))
            new{T,X}(c, Ref(m′),  Ref(n))
        end
    end
    # non copying version assumes trimmed coeffs
    function LaurentPolynomial{T,X}(::Val{false}, coeffs::AbstractVector{T},
                                    m::Integer=0) where {T, X}
        new{T,X}(coeffs, Ref(m), Ref(m + length(coeffs) - 1))
    end
end

@register LaurentPolynomial

## constructors
function LaurentPolynomial{T}(coeffs::AbstractVector{S}, m::Int, var::SymbolLike=Var(:x)) where {
    T, S <: Number}
    LaurentPolynomial{T,Symbol(var)}(T.(coeffs), m)
end

function LaurentPolynomial{T}(coeffs::AbstractVector{T}, var::SymbolLike=Var(:x)) where {T}
    LaurentPolynomial{T, Symbol(var)}(coeffs, 0)
end

function LaurentPolynomial(coeffs::AbstractVector{T}, m::Int, var::SymbolLike=Var(:x)) where {T}
    LaurentPolynomial{T, Symbol(var)}(coeffs, m)
end




##
## conversion
##

# LaurentPolynomial is a wider collection than other standard basis polynomials.
Base.promote_rule(::Type{P},::Type{Q}) where {T, X, P <: LaurentPolynomial{T,X}, S, Q <: StandardBasisPolynomial{S, X}} = LaurentPolynomial{promote_type(T, S), X}


Base.promote_rule(::Type{Q},::Type{P}) where {T, X, P <: LaurentPolynomial{T,X}, S, Q <: StandardBasisPolynomial{S,X}} =
    LaurentPolynomial{promote_type(T, S),X}

# need to add p.m[], so abstract.jl method isn't sufficent
# XXX unlike abstract.jl, this uses Y variable in conversion; no error
# Used in DSP.jl
function Base.convert(::Type{LaurentPolynomial{S,Y}}, p::LaurentPolynomial{T,X}) where {T,X,S,Y}
    LaurentPolynomial{S,Y}(p.coeffs, p.m[])
end

# work around for non-applicable convert(::Type{<:P}, p::P{T,X}) in abstract.jl
struct OffsetCoeffs{V}
    coeffs::V
    m::Int
end

_coeffs(p::LaurentPolynomial) = OffsetCoeffs(p.coeffs, p.m[])
function LaurentPolynomial{T,X}(p::OffsetCoeffs) where {T, X}
    LaurentPolynomial{T,X}(p.coeffs, p.m)
end

function Base.convert(::Type{P}, q::StandardBasisPolynomial{S}) where {P <:LaurentPolynomial,S}

     T = _eltype(P, q)
     X = indeterminate(P, q)
     ⟒(P){T,X}([q[i] for i in eachindex(q)], firstindex(q))
 end

function Base.convert(::Type{P}, q::AbstractPolynomial{T,X}) where {T,X,P <:LaurentPolynomial}
    convert(P, convert(Polynomial, q))
end

# Ambiguity, issue #435
function Base.convert(𝑷::Type{P}, p::ArnoldiFit{T, M, X}) where {P<:LaurentPolynomial, T, M, X}
    convert(𝑷, convert(Polynomial, p))
end

##
## generic functions
##

function Base.inv(p::LaurentPolynomial{T, X}) where {T, X}
    m,n =  (extrema∘degreerange)(p)
    m != n && throw(ArgumentError("Only monomials can be inverted"))
    cs = [1/p for p in p.coeffs]
    LaurentPolynomial{eltype(cs), X}(cs, -m)
end

Base.numerator(p::LaurentPolynomial) = numerator(convert(RationalFunction, p))
Base.denominator(p::LaurentPolynomial) = denominator(convert(RationalFunction, p))

##
## changes to common.jl mostly as the range in the type is different
##
Base.:(==)(p1::LaurentPolynomial, p2::LaurentPolynomial) =
    check_same_variable(p1, p2) && (degreerange(chop!(p1)) == degreerange(chop!(p2))) && (coeffs(p1) == coeffs(p2))
Base.hash(p::LaurentPolynomial, h::UInt) = hash(indeterminate(p), hash(degreerange(p), hash(coeffs(p), h)))

isconstant(p::LaurentPolynomial) = iszero(lastindex(p)) && iszero(firstindex(p))

basis(P::Type{<:LaurentPolynomial{T,X}}, n::Int) where {T,X} = LaurentPolynomial{T,X}(ones(T,1), n)

# like that in common, only return zero if idx < firstindex(p)
function Base.getindex(p::LaurentPolynomial{T}, idx::Int) where {T}
    m,M = firstindex(p), lastindex(p)
    m <= idx <= M || return zero(T)
    p.coeffs[idx-m+1]
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

minimumexponent(::Type{<:LaurentPolynomial}) = typemin(Int)
minimumexponent(p::LaurentPolynomial) = p.m[]
Base.firstindex(p::LaurentPolynomial) = minimumexponent(p)
degree(p::LaurentPolynomial) = p.n[]


_convert(p::P, as) where {T,X,P <: LaurentPolynomial{T,X}} = ⟒(P)(as, firstindex(p), Var(X))

## chop!
# trim  from *both* ends
function chop!(p::LaurentPolynomial{T};
               rtol::Real = Base.rtoldefault(real(T)),
               atol::Real = 0,) where {T}

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

    cs = copy(coeffs(p))
    rng = m-m0+1:n-m0+1
    resize!(p.coeffs, length(rng))
    p.coeffs[:] = cs[rng]
    isempty(p.coeffs) && push!(p.coeffs,zero(T))
    p.m[], p.n[] = m, max(m,n)

    p

end

# use unicode exponents. XXX modify printexponent to always use these?
function showterm(io::IO, ::Type{<:LaurentPolynomial}, pj::T, var, j, first::Bool, mimetype) where {T}
    if iszero(pj) return false end
    pj = printsign(io, pj, first, mimetype)
    if !(hasone(T) && pj == one(T) && !(showone(T) || j == 0))
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

```jldoctest laurent
julia> using Polynomials;

julia> z = variable(LaurentPolynomial, :z)
LaurentPolynomial(1.0*z)

julia> p = LaurentPolynomial([im, 1+im, 2 + im], -1, :z)
LaurentPolynomial(im*z⁻¹ + (1 + im) + (2 + im)z)

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

```jldoctest laurent
julia> using Polynomials;

julia> z = variable(LaurentPolynomial, :z)
LaurentPolynomial(z)

julia> h = LaurentPolynomial([1,1], -1, :z)
LaurentPolynomial(z⁻¹ + 1)

julia> Polynomials.paraconj(h)(z) ≈ 1 + z ≈ LaurentPolynomial([1,1], 0, :z)
true

julia> h = LaurentPolynomial([3,2im,1], -2, :z)
LaurentPolynomial(3*z⁻² + 2im*z⁻¹ + 1)

julia> Polynomials.paraconj(h)(z) ≈ 1 - 2im*z + 3z^2 ≈ LaurentPolynomial([1, -2im, 3], 0, :z)
true

julia> Polynomials.paraconj(h)(z) ≈ (conj ∘ h ∘ conj ∘ inv)(z)
true
"""
function paraconj(p::LaurentPolynomial)
    cs = p.coeffs
    ds = adjoint.(cs)
    n = degree(p)
    LaurentPolynomial(reverse(ds), -n, indeterminate(p))
end

"""
    cconj(p)

Conjugation of a polynomial with respect to the imaginary axis.

The `cconj` of a polynomial, `p̃`, conjugates the coefficients and applies `s -> -s`. That is `cconj(p)(s) = conj(p)(-s)`.

This satisfies for *imaginary* `s`: `conj(p(s)) = p̃(s) = (conj ∘ p)(s) = cconj(p)(s) `

[ref](https://github.com/hurak/PolynomialEquations.jl#symmetrix-conjugate-equation-continuous-time-case)

Examples:
```jldoctest laurent
julia> using Polynomials;

julia> s = 2im
0 + 2im

julia> p = LaurentPolynomial([im,-1, -im, 1], 1, :s)
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
    LaurentPolynomial(ps, m, indeterminate(p))
end



##
## ----
##


# evaluation uses `evalpoly`
function evalpoly(x::S, p::LaurentPolynomial{T}) where {T,S}
    xᵐ = firstindex(p) < 0 ? inv(x)^(abs(firstindex(p))) : (x/1)^firstindex(p) # make type stable
    return EvalPoly.evalpoly(x, p.coeffs) * xᵐ
end



# scalar operations
# needed as standard-basis defn. assumes basis 1, x, x², ...
function Base.:+(p::LaurentPolynomial{T,X}, c::S) where {T, X, S <: Number}
    R = promote_type(T,S)
    q = LaurentPolynomial{R,X}(p.coeffs, firstindex(p))
    q[0] += c
    q
end

##
## Poly +, - and  *
## uses some ideas from https://github.com/jmichel7/LaurentPolynomials.jl/blob/main/src/LaurentPolynomials.jl for speedups
Base.:+(p1::LaurentPolynomial, p2::LaurentPolynomial) = add_sub(+, p1, p2)
Base.:-(p1::LaurentPolynomial, p2::LaurentPolynomial) = add_sub(-, p1, p2)

function add_sub(op, p1::P, p2::Q) where {T, X, P <: LaurentPolynomial{T,X},
                                          S, Y, Q <: LaurentPolynomial{S,Y}}

    isconstant(p1) && return op(constantterm(p1), p2)
    isconstant(p2) && return op(p1, constantterm(p2))
    assert_same_variable(X, Y)

    m1,n1 = (extrema ∘ degreerange)(p1)
    m2,n2 = (extrema ∘ degreerange)(p2)
    m, n = min(m1,m2), max(n1, n2)

    R = promote_type(T,S)
    as = zeros(R, length(m:n))

    d = m1 - m2
    d1, d2 = m1 > m2 ? (d,0) : (0, -d)

    for (i, pᵢ) ∈ pairs(p1.coeffs)
        @inbounds as[d1 + i] = pᵢ
    end
    for (i, pᵢ) ∈ pairs(p2.coeffs)
        @inbounds as[d2 + i] = op(as[d2+i], pᵢ)
    end

    m = _laurent_chop!(as, m)
    isempty(as) && return zero(LaurentPolynomial{R,X})
    q = LaurentPolynomial{R,X}(Val(false), as, m)
    return q

end

function Base.:*(p1::P, p2::Q) where {T,X,P<:LaurentPolynomial{T,X},
                                      S,Y,Q<:LaurentPolynomial{S,Y}}

    isconstant(p1) && return constantterm(p1) * p2
    isconstant(p2) && return p1 * constantterm(p2)
    assert_same_variable(X, Y)

    m1,n1 = (extrema ∘ degreerange)(p1)
    m2,n2 = (extrema ∘ degreerange)(p2)
    m,n = m1 + m2, n1+n2

    R = promote_type(T,S)
    as = zeros(R, length(m:n))

    for (i, p₁ᵢ) ∈ pairs(p1.coeffs)
        for (j, p₂ⱼ) ∈ pairs(p2.coeffs)
            @inbounds as[i+j-1] += p₁ᵢ * p₂ⱼ
        end
    end

    m = _laurent_chop!(as, m)

    isempty(as) && return zero(LaurentPolynomial{R,X})
    p = LaurentPolynomial{R,X}(Val(false), as, m)

    return p
end

function _laurent_chop!(as, m)
    while !isempty(as)
        if iszero(first(as))
            m += 1
            popfirst!(as)
        else
            break
        end
    end
    while !isempty(as)
        if iszero(last(as))
            pop!(as)
        else
            break
        end
    end
    m
end

function scalar_mult(p::LaurentPolynomial{T,X}, c::Number) where {T,X}
    LaurentPolynomial(p.coeffs .* c, p.m[], Var(X))
end
function scalar_mult(c::Number, p::LaurentPolynomial{T,X}) where {T,X}
    LaurentPolynomial(c .* p.coeffs, p.m[], Var(X))
end


function integrate(p::P) where {T, X, P <: LaurentBasisPolynomial{T, X}}

    R = typeof(constantterm(p) / 1)
    Q = ⟒(P){R,X}

    hasnan(p) && return Q([NaN])
    iszero(p) && return zero(Q)

    ∫p = zero(Q)
    for (k, pₖ) ∈ pairs(p)
        iszero(pₖ) && continue
        k == -1 && throw(ArgumentError("Can't integrate Laurent polynomial with  `x⁻¹` term"))
        ∫p[k+1] = pₖ/(k+1)
    end
    ∫p
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
julia> using Polynomials;

julia> p = LaurentPolynomial([24,10,-15,0,1],-2,:z)
LaurentPolynomial(24*z⁻² + 10*z⁻¹ - 15 + z²)

julia> roots(p)
4-element Vector{Float64}:
 -3.999999999999999
 -0.9999999999999994
  1.9999999999999998
  2.9999999999999982
```
"""
function  roots(p::P; kwargs...)  where  {T, X, P <: LaurentPolynomial{T, X}}
    c = coeffs(p)
    r = degreerange(p)
    d = r[end] - min(0, r[1]) + 1    # Length of the coefficient vector, taking into consideration
                                     # the case when the lower degree is strictly positive
                                     # (like p=3z^2).
    z = zeros(T, d)                  # Reserves space for the coefficient vector.
    z[max(0, r[1]) + 1:end] = c      # Leaves the coeffs of the lower powers as zeros.
    a = Polynomial{T,X}(z)           # The root is then the root of the numerator polynomial.
    return roots(a; kwargs...)
end

##
## d/dx, ∫
##
function derivative(p::P, order::Integer = 1) where {T, X, P<:LaurentPolynomial{T,X}}

    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return p

    hasnan(p) && return ⟒(P)(T[NaN], 0, X)

    m,n = (extrema ∘ degreerange)(p)
    m = m - order
    n = n - order
    as =  zeros(T, length(m:n))

    for (k, pₖ) in pairs(p)
        iszero(pₖ) && continue
        idx = 1 + k - order - m
        if 0 ≤ k ≤ order - 1
            as[idx] = zero(T)
        else
            as[idx] = reduce(*, (k - order + 1):k, init = pₖ)
        end
    end

    chop!(LaurentPolynomial{T,X}(as, m))

end


function Base.gcd(p::LaurentPolynomial{T,X}, q::LaurentPolynomial{T,Y}, args...; kwargs...) where {T,X,Y}
    mp, Mp = (extrema ∘ degreerange)(p)
    mq, Mq = (extrema ∘ degreerange)(q)
    if mp < 0 || mq < 0
        throw(ArgumentError("GCD is not defined when there are `x⁻ⁿ` terms"))
    end

    degree(p) == 0 && return iszero(p) ? q : one(q)
    degree(q) == 0 && return iszero(q) ? p : one(p)
    assert_same_variable(p,q)

    pp, qq = convert(Polynomial, p), convert(Polynomial, q)
    u = gcd(pp, qq, args..., kwargs...)
    return LaurentPolynomial(coeffs(u), X)
end
