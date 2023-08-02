# Dense + StandardBasis

"""
    LaurentPolynomial{T,X}(coeffs::AbstractVector, [m::Integer = 0], [var = :x])

A [Laurent](https://en.wikipedia.org/wiki/Laurent_polynomial) polynomial is of the form `a_{m}x^m + ... + a_{n}x^n` where `m,n` are  integers (not necessarily positive) with ` m <= n`.

The `coeffs` specify `a_{m}, a_{m-1}, ..., a_{n}`.
The argument `m` represents the lowest exponent of the variable in the series, and is taken to be zero by default.

Laurent polynomials and standard basis polynomials promote to  Laurent polynomials. Laurent polynomials may be  converted to a standard basis  polynomial when `m >= 0`
.

Integration will fail if there is a `x⁻¹` term in the polynomial.

!!! note
    `LaurentPolynomial` is an alias for `MutableDensePolynomial{StandardBasis}`.

!!! note
    `LaurentPolynomial` is axis-aware, unlike the other polynomial types in this package.

# Examples:
```jldoctest laurent
julia> using Polynomials

julia> P = LaurentPolynomial;

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
const LaurentPolynomial = MutableDenseLaurentPolynomial{StandardBasis}
export LaurentPolynomial

_typealias(::Type{P}) where {P<:LaurentPolynomial} = "LaurentPolynomial"

# how to show term. Only needed here to get unicode exponents to match old LaurentPolynomial.jl type
function showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {T, P<:MutableDenseLaurentPolynomial{StandardBasis,T}}
    if _iszero(pj) return false end

    pj = printsign(io, pj, first, mimetype)
    if hasone(T)
        if !(_isone(pj) && !(showone(T) || j == 0))
            printcoefficient(io, pj, j, mimetype)
        end
    else
        printcoefficient(io, pj, j, mimetype)
    end

    iszero(j) && return true
    printproductsign(io, pj, j, mimetype)
    print(io, indeterminate(P))
    j == 1 && return true
    unicode_exponent(io, j) # print(io, exponent_text(j, mimetype))
    return true
end


function evalpoly(c, p::LaurentPolynomial{T,X}) where {T,X}
    iszero(p) && return zero(T) * zero(c)
    EvalPoly.evalpoly(c, p.coeffs) * c^p.order[]
end

# scalar add
function scalar_add(c::S, p::LaurentPolynomial{T,X}) where {S, T, X}
    R = promote_type(T,S)
    P = LaurentPolynomial{R,X}

    iszero(p) && return P([c], 0)
    iszero(c) && return convert(P, p)

    a,b = firstindex(p), lastindex(p)
    a′ = min(0, a)
    b′ = max(0, b)
    cs = _zeros(p, zero(first(p.coeffs) + c), length(a′:b′))
    o = offset(p) + a - a′
    for (i, cᵢ) ∈ pairs(p)
        cs[i + o] = cᵢ
    end
    cs[0 + o] += c
    iszero(last(cs)) && (cs = trim_trailing_zeros(cs))
    P(Val(false), cs, a′)
end

function ⊗(p::LaurentPolynomial{T,X},
           q::LaurentPolynomial{S,X}) where {T,S,X}
    # simple convolution
    # This is ⊗(P,p,q) from polynomial standard-basis
    R = promote_type(T,S)
    P = LaurentPolynomial{R,X}

    iszero(p) && return zero(P)
    iszero(q) && return zero(P)

    a₁, a₂ = firstindex(p), firstindex(q)
    b₁, b₂ = lastindex(p), lastindex(q)
    a, b = a₁ + a₂, b₁ + b₂

    z = zero(first(p) * first(q))
    cs = _zeros(p, z, length(a:b))

    # convolve and shift order
    @inbounds for (i, pᵢ) ∈ enumerate(p.coeffs)
        for (j, qⱼ) ∈ enumerate(q.coeffs)
            ind = i + j - 1
            cs[ind] = muladd(pᵢ, qⱼ, cs[ind])
        end
    end
    if iszero(last(cs))
        cs = trim_trailing_zeros(cs)
    end
    P(Val(false), cs, a)
end

function derivative(p::LaurentPolynomial{T,X}) where {T,X}

    N = lastindex(p) - firstindex(p) + 1
    R = promote_type(T, Int)
    P = LaurentPolynomial{R,X}
    hasnan(p) && return P(zero(T)/zero(T)) # NaN{T}
    iszero(p) && return P(0*p[0])

    ps = p.coeffs
    cs = [i*pᵢ for (i,pᵢ) ∈ pairs(p)]
    return P(cs, p.order[]-1)
end

# LaurentPolynomials have `inv` defined for monomials
function Base.inv(p::LaurentPolynomial)
    m,n =  firstindex(p), lastindex(p)
    m != n && throw(ArgumentError("Only monomials can be inverted"))
    cs = [1/p for p in p.coeffs]
    LaurentPolynomial{eltype(cs), indeterminate(p)}(cs, -m)
end

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
function paraconj(p::LaurentPolynomial{T,X}) where {T,X}
    cs = p.coeffs
    ds = adjoint.(cs)
    n = degree(p)
    LaurentPolynomial{T,X}(reverse(ds), -n)
end
paraconj(::AbstractPolynomial) = throw(ArgumentError("`paraconj` not defined for this polynomial type"))

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
function cconj(p::LaurentPolynomial{T,X}) where {T,X}
    ps = conj.(coeffs(p))
    m,n = firstindex(p), lastindex(p)
    for i in m:n
        if isodd(i)
            ps[i+1-m] *= -1
        end
    end
    LaurentPolynomial{T,X}(ps, m)
end
cconj(::AbstractPolynomial) = throw(ArgumentError("`cconj` not defined for this polynomial type"))


function  roots(p::P; kwargs...)  where  {T, X, P <:LaurentPolynomial{T,X}}
    return roots(convert(Polynomial, numerator(p)), kwargs...)
end
