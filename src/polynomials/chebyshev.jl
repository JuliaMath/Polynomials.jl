# Example of using mutable dense container with a different basis
struct ChebyshevTBasis <: AbstractBasis end

# This is the same as ChebyshevT
#const ChebyshevT = MutableDensePolynomial{ChebyshevTBasis}
#export ChebyshevT

basis_symbol(::Type{<:AbstractUnivariatePolynomial{ChebyshevTBasis}}) = "T"

function Base.convert(P::Type{<:Polynomial}, ch::MutableDensePolynomial{ChebyshevTBasis})

    T = _eltype(P,ch)
    X = indeterminate(P,ch)
    Q = ⟒(P){T,X}

    d = lastindex(ch)
    if d ≤ 1 # T₀, T₁ = 1, x
        return Q(coeffs(ch))
    end

    c0 = Q(ch[end - 1])
    c1 = Q(ch[end])
    x = variable(Q)
    @inbounds for i in d:-1:2
        tmp = c0
        c0 = Q(ch[i - 2]) - c1
        c1 = tmp + c1 * x * 2
    end
    return c0 + c1 * x
end

function Base.convert(C::Type{<:MutableDensePolynomial{ChebyshevTBasis}}, p::Polynomial)
    x = variable(C)
    isconstant(p) || assert_same_variable(indeterminate(x),indeterminate(p))
    p(x)
end

# lowest degree is always 0
laurenttype(::Type{<:MutableDensePolynomial{ChebyshevTBasis}}) = Val(false)
minimumexponent(::Type{<:MutableDensePolynomial{ChebyshevTBasis}}) = 0
domain(::Type{<:MutableDensePolynomial{ChebyshevTBasis}}) = Interval(-1, 1)

constantterm(p::MutableDensePolynomial{ChebyshevTBasis}) = p(0)
function Base.one(::Type{P}) where {P<:MutableDensePolynomial{ChebyshevTBasis}}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}(ones(T,1))
end
function variable(::Type{P}) where {P<:MutableDensePolynomial{ChebyshevTBasis}}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}([zero(T), one(T)])
end

"""
    (::MutableDensePolynomial{ChebyshevTBasis})(x)

Evaluate the Chebyshev polynomial at `x`. If `x` is outside of the domain of [-1, 1], an error will be thrown. The evaluation uses Clenshaw Recursion.

# Examples
```jldoctest ChebyshevT
julia> using Polynomials

julia> c = ChebyshevT([2.5, 1.5, 1.0])
ChebyshevT(2.5⋅T_0(x) + 1.5⋅T_1(x) + 1.0⋅T_2(x))

julia> c(0)
1.5

julia> c.(-1:0.5:1)
5-element Vector{Float64}:
 2.0
 1.25
 1.5
 2.75
 5.0
```
"""
function evalpoly(x::S, ch::MutableDensePolynomial{ChebyshevTBasis}) where {S}
    x ∉ domain(ch) && throw(ArgumentError("$x outside of domain"))
    evalpoly(x, ch, false)
end

function evalpoly(x::AbstractPolynomial, ch::MutableDensePolynomial{ChebyshevTBasis})
    evalpoly(x, ch, false)
end

# no checking, so can be called directly through any third argument
function evalpoly(x::S, ch::MutableDensePolynomial{ChebyshevTBasis,T}, checked) where {T,S}
    R = promote_type(T, S)
    length(ch) == 0 && return zero(R)
    length(ch) == 1 && return R(ch[0])
    c0 = ch[end - 1]
    c1 = ch[end]
    @inbounds for i in lastindex(ch) - 2:-1:0
        c0, c1 = ch[i] - c1, c0 + c1 * 2x
    end
    return R(c0 + c1 * x)
end

# scalar +
function scalar_add(c::S, p::MutableDensePolynomial{B,T,X}) where {B<:ChebyshevTBasis,T,X, S<:Scalar}
    R = promote_type(T,S)
    cs = collect(R, values(p))
    cs[1] += c
    MutableDensePolynomial{ChebyshevTBasis,R,X}(cs)
end

# product
function ⊗(p1::MutableDensePolynomial{B,T,X}, p2::MutableDensePolynomial{B,T,X}) where {B<:ChebyshevTBasis,T,X}
    z1 = _c_to_z(coeffs(p1))
    z2 = _c_to_z(coeffs(p2))
    prod = fastconv(z1, z2)
    cs = _z_to_c(prod)
    ret = MutableDensePolynomial{ChebyshevTBasis}(cs,X)
    return ret
end

function derivative(p::P, order::Integer = 1) where {B<:ChebyshevTBasis,T,X,P<:MutableDensePolynomial{B,T,X}}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    R  = eltype(one(T)/1)
    Q = MutableDensePolynomial{ChebyshevTBasis,R,X}
    order == 0 && return convert(Q, p)
    hasnan(p) && return Q(R[NaN])
    order > length(p) && return zero(Q)


    q =  convert(P{R,X}, copy(p))
    n = length(p)
    der = Vector{R}(undef, n)

    for j in n:-1:3
        der[j] = 2j * q[j]
        q[j - 2] += j * q[j] / (j - 2)
    end
    if n > 1
        der[2] = 4q[2]
    end
    der[1] = q[1]

    pp = Q(der)
    return order > 1 ?  derivative(pp, order - 1) : pp

end

function integrate(p::P) where {B<:ChebyshevTBasis,T,X,P<:MutableDensePolynomial{B,T,X}}
    R = eltype(one(T) / 1)
    Q = MutableDensePolynomial{B,R,X}
    if hasnan(p)
        return Q([NaN])
    end
    n = length(p)
    if n == 1
        return Q([zero(R), p[0]])
    end
    a2 = Vector{R}(undef, n + 1)
    a2[1] = zero(R)
    a2[2] = p[0]
    a2[3] = p[1] / 4
    @inbounds for i in 2:n - 1
        a2[i + 2] = p[i] / (2 * (i + 1))
        a2[i] -= p[i] / (2 * (i - 1))
    end

    return Q(a2)
end

function vander(P::Type{<:MutableDensePolynomial{ChebyshevTBasis}}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    if n > 0
        A[:, 2] .= x
        @inbounds for i in 3:n + 1
            A[:, i] .= A[:, i - 1] .* 2x .- A[:, i - 2]
        end
    end
    return A
end

function companion(p::MutableDensePolynomial{ChebyshevTBasis,T}) where T
    d = length(p) - 1
    d < 1 && throw(ArgumentError("Series must have degree greater than 1"))
    d == 1 && return diagm(0 => [-p[0] / p[1]])
    R = eltype(one(T) / one(T))

    scl = vcat(1.0, fill(R(√0.5), d - 1))

    diag = vcat(√0.5, fill(R(0.5), d - 2))
    comp = diagm(1 => diag,
                 -1 => diag)
    monics = coeffs(ps) ./ coeffs(p)[end]
    comp[:, end] .-= monics[1:d] .* scl ./ scl[end] ./ 2
    return R.(comp)
end

function Base.divrem(num::MutableDensePolynomial{ChebyshevTBasis}{T,X},
                     den::MutableDensePolynomial{ChebyshevTBasis}{S,Y}) where {T,X,S,Y}
    assert_same_variable(num, den)
    n = length(num) - 1
    m = length(den) - 1

    R = typeof(one(T) / one(S))
    P = MutableDensePolynomial{ChebyshevTBasis}{R,X}

    if n < m
        return zero(P), convert(P, num)
    elseif m == 0
        den[0] ≈ 0 && throw(DivideError())
        return num ./ den[end], zero(P)
    end

    znum = _c_to_z(coeffs(num))
    zden = _c_to_z(coeffs(den))
    quo, rem = _z_division(znum, zden)
    q_coeff = _z_to_c(quo)
    r_coeff = _z_to_c(rem)
    return P(q_coeff), P(r_coeff)
end
