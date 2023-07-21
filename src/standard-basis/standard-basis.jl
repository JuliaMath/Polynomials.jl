struct StandardBasis <: AbstractBasis end

function print_basis(io::IO, p::AbstractUnivariatePolynomial{<:StandardBasis}, i)
    print(io, X)
    print_unicode_exponent(io, i)
end

# XXX For now need 3 convert methods for standard basis
function Base.convert(P::Type{PP}, q::Q) where {B<:StandardBasis, PP <: AbstractUnivariatePolynomial{B}, Q<:AbstractUnivariatePolynomial{B}}
    if isa(q, PP)
        return q
    else
        minimumexponent(P) <= minimumexponent(q) ||
            throw(ArgumentError("a $P can not have a minimum exponent of $(minimumexponent(q))"))
        T = _eltype(P,q)
        X = indeterminate(P,q)
        return ⟒(P){T,X}([q[i] for i in eachindex(q)], firstindex(q))
    end
end
function Base.convert(P::Type{PP}, q::Q) where {PP <: StandardBasisPolynomial, B<:StandardBasis,T,X, Q<:AbstractUnivariatePolynomial{B,T,X}}
    minimumexponent(P) > firstindex(q) &&
        throw(ArgumentError("Degree of polynomial less than minimum degree of polynomial type $(⟒(P))"))

    isa(q, PP) && return p
    T′ = _eltype(P,q)
    X′ = indeterminate(P,q)
    if firstindex(q) >= 0
        cs = [q[i] for i ∈ 0:lastindex(q)]
        o = 0
    else
        cs = [q[i] for i ∈ eachindex(q)]
        o = firstindex(q)
    end
    ⟒(P){T′,X′}(cs, o)
end
function Base.convert(P::Type{PP}, q::Q) where {B<:StandardBasis, PP <: AbstractUnivariatePolynomial{B}, Q<:StandardBasisPolynomial}
    isa(q, PP) && return p
    T = _eltype(P,q)
    X = indeterminate(P,q)
    ⟒(P){T,X}([q[i] for i in eachindex(q)], firstindex(q))
end

Base.one(p::P) where {B<:StandardBasis,T,X, P <: AbstractUnivariatePolynomial{B,T,X}} = ⟒(P){T,X}([one(p[0])])
Base.one(::Type{P}) where {B<:StandardBasis,T,P <: AbstractUnivariatePolynomial{B,T}} = ⟒(P){T}(ones(T,1), indeterminate(P))
Base.one(::Type{P}) where {B<:StandardBasis,T,X, P <: AbstractUnivariatePolynomial{B,T,X}} = ⟒(P){T,X}(ones(T,1))

variable(P::Type{<:AbstractUnivariatePolynomial{B,T,X}}) where {B<:StandardBasis,T,X} =  basis(P, 1)

basis(P::Type{<:AbstractUnivariatePolynomial{B, T, X}}, i::Int) where {B<:StandardBasis,T,X} = P(ones(T,1), i)

constantterm(p::AbstractUnivariatePolynomial{B}) where {B <: StandardBasis} = p[0]

domain(::Type{P}) where {B <: StandardBasis, P <: AbstractUnivariatePolynomial{B}} = Interval{Open,Open}(-Inf, Inf)

mapdomain(::Type{P}, x::AbstractArray) where  {B <: StandardBasis, P <: AbstractUnivariatePolynomial{B}} = x

## Multiplication
# special cases are faster
function ⊗(p::AbstractUnivariatePolynomial{B,T,X},
           q::AbstractUnivariatePolynomial{B,S,X}) where {B <: StandardBasis, T,S,X}
    # simple convolution with order shifted
    XXX()
end

function integrate(p::AbstractUnivariatePolynomial{B,T,X}) where {B <: StandardBasis,T,X}

    iszero(p) && return p/1
    N = lastindex(p) - firstindex(p) + 1
    R = typeof(one(T)/1)
    z = zero(R)
    P = ⟒(p){R,X}
    hasnan(p) && return  P(zero(T)/zero(T)) # NaN{T}

    cs = _zeros(p, z, N+1)
    os =  offset(p)
    @inbounds for (i, cᵢ) ∈ pairs(p)
        i == -1 && (iszero(cᵢ) ? continue : throw(ArgumentError("Laurent polynomial with 1/x term")))
        #cs[i + os] = cᵢ / (i+1)
        cs = _set(cs, i + 1 + os,  cᵢ / (i+1))
    end
    P(cs, firstindex(p))
end

function Base.divrem(num::P, den::Q) where {B<:StandardBasis,
                                            T, P <: AbstractUnivariatePolynomial{B,T},
                                            S, Q <: AbstractUnivariatePolynomial{B,S}}

    assert_same_variable(num, den)
    @assert ⟒(P) == ⟒(Q)

    X = indeterminate(num)
    R = Base.promote_op(/, T, S)
    PP = ⟒(P){R,X}


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
    return PP(q_coeff), PP(r_coeff)

end

## XXX This needs resolving!
## XXX copy or pass along to other system for now where things are defined fro StandardBasisPolynomial
function vander(p::Type{<:P}, x::AbstractVector{T}, degs) where {B<:StandardBasis, P<:AbstractUnivariatePolynomial{B}, T <: Number}
    vander(StandardBasisPolynomial, x, degs)
end

function LinearAlgebra.cond(p::P, x) where {B<:StandardBasis, P<:AbstractUnivariatePolynomial{B}}
    p̃ = map(abs, p)
    p̃(abs(x))/ abs(p(x))
end

function ngcd(p::P, q::Q,
              args...;
              kwargs...) where {B <: StandardBasis,
                                T,X,P<:AbstractUnivariatePolynomial{B,T,X},
                                S,Y,Q<:AbstractUnivariatePolynomial{B,S,Y}}
    ngcd(PnPolynomial(p.coeffs), PnPolynomial(q.coeffs), args...; kwargs...)
end

# XXX p.coeffs isn't right
function Multroot.multroot(p::AbstractUnivariatePolynomial{B}, args...;
                           kwargs...) where {B<:StandardBasis}
    cs = coeffs(p)
    if firstindex(p) > 0
        cs = vcat(zeros(firstindex(p)), cs)
    elseif firstindex(p) < 0
        @warn "Laurent Polynomial; finding values after factoring out leading term"
    end
    Multroot.multroot(Polynomial(cs), args...; kwargs...)
end

Polynomials.Multroot.pejorative_root(q::AbstractUnivariatePolynomial{<:StandardBasis}, zs::Vector{S}, ls; kwargs...) where {S} =
    Polynomials.Multroot.pejorative_root(convert(Polynomial, q), zs, ls; kwargs...)

Polynomials.Multroot.stats(q::AbstractUnivariatePolynomial{<:StandardBasis}, zs::Vector{S}, ls; kwargs...) where {S} =
    Polynomials.Multroot.stats(convert(Polynomial, q), zs, ls; kwargs...)

function fit(::Type{P},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::Int;
             kwargs...) where {T, P<:AbstractUnivariatePolynomial{<:StandardBasis}}
    convert(P, fit(Polynomial, x, y, deg; kwargs...))
end

# for one ambigous test!
function fit(::Type{P},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::AbstractVector;
             kwargs...) where {T, P<:AbstractUnivariatePolynomial{<:StandardBasis}}
    convert(P, fit(Polynomial, x, y, deg; kwargs...))
end

function fit(::Type{P},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::AbstractVector,
             cs::Dict;
             kwargs...) where {T, P<:AbstractUnivariatePolynomial{<:StandardBasis}}
    convert(P, fit(Polynomial, x, y, deg, cs; kwargs...))
end

# new constructors taking order in second position
#SparsePolynomial{T,X}(coeffs::AbstractVector{S}, ::Int) where {T, X, S} = SparsePolynomial{T,X}(coeffs)
Polynomial{T, X}(coeffs::AbstractVector{S},order::Int) where {T, X, S} = Polynomial{T,X}(coeffs)
#ImmutablePolynomial{T,X}(coeffs::AbstractVector{S}, ::Int)  where {T,X,S} = ImmutablePolynomial{T,X}(coeffs)
FactoredPolynomial{T,X}(coeffs::AbstractVector{S}, order::Int) where {T,S,X} = FactoredPolynomial{T,X}(coeffs)
PnPolynomial{T, X}(coeffs::AbstractVector, order::Int) where {T, X} = PnPolynomial(coeffs) # for generic programming
