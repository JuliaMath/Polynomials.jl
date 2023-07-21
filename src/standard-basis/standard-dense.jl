# Dense + StandardBasis


function evalpoly(c, p::MutableDensePolynomial{StandardBasis,T,X}) where {T,X}
    iszero(p) && return zero(T) * zero(c)
    EvalPoly.evalpoly(c, p.coeffs) * c^p.order
end


# scalar add
function scalar_add(c::S, p:: MutableDensePolynomial{StandardBasis,T,X}) where {S, T, X}
    R = promote_type(T,S)
    P =  MutableDensePolynomial{StandardBasis,R,X}

    iszero(p) && return P([c], 0)
    iszero(c) && return convert(P, p)

    a,b = firstindex(p), lastindex(p)
    a′ = min(0,a)
    cs = _zeros(p, zero(first(p.coeffs)+c), length(a′:b))
    o = offset(p) + a - a′
    for (i, cᵢ) ∈ pairs(p)
        cs[i+o] = cᵢ
    end
    cs[0+o] += c
    iszero(last(cs)) && (cs = trim_trailing_zeros(cs))
    P(Val(false), cs, a′)
end



function ⊗(p:: MutableDensePolynomial{StandardBasis,T,X},
           q:: MutableDensePolynomial{StandardBasis,S,X}) where {T,S,X}
    # simple convolution
    # This is ⊗(P,p,q) from polynomial standard-basis
    R = promote_type(T,S)
    P =  MutableDensePolynomial{StandardBasis,R,X}

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


function derivative(p::MutableDensePolynomial{B,T,X}) where {B<:StandardBasis,T,X}

    N = lastindex(p) - firstindex(p) + 1
    R = promote_type(T, Int)
    P = ⟒(p){R,X}
    hasnan(p) && return P(zero(T)/zero(T)) # NaN{T}
    iszero(p) && return P(0*p[0])

    ps = p.coeffs
    cs = [i*pᵢ for (i,pᵢ) ∈ pairs(p)]
    return P(cs, p.order-1)
end

# LaurentPolynomials have `inv` defined for monomials
function Base.inv(p::MutableDensePolynomial{StandardBasis})
    m,n =  firstindex(p), lastindex(p)
    m != n && throw(ArgumentError("Only monomials can be inverted"))
    cs = [1/p for p in p.coeffs]
    LaurentPolynomial{eltype(cs), indeterminate(p)}(cs, -m)
end

## XXX ----
#const Polynomial = MutableDensePolynomial{StandardBasis}
#export Polynomial

const LaurentPolynomial = MutableDensePolynomial{StandardBasis}
export LaurentPolynomial

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



# resolve ambiguity
# function Base.convert(::Type{P}, q::Q) where {T, X, P<:LaurentPolynomial, B<:StandardBasis, Q<:AbstractUnivariatePolynomial{B, T, X}}
#     p = convert(Polynomial, q)
#     LaurentPolynomial{eltype(p), indeterminate(p)}(p.coeffs, 0)
# end

## ----
## XXX needs to be incorporated if Polynomial =  MutableDensePolynomial{StandardBasis}
function  roots(p::P; kwargs...)  where  {T, X, P <: MutableDensePolynomial{StandardBasis,T,X}}
    iszero(p) && return float(T)[]
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
