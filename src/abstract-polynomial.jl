"""
    Abstract type for polynomials with an explicit basis.
"""
abstract type AbstractUnivariatePolynomial{B, T, X} <: AbstractPolynomial{T,X} end
abstract type AbstractBasis end

## idea is vector space stuff (scalar_add, scalar_mult, vector +/-, ^) goes here
## connection (convert, transform) is specific to a basis (storage)
## ⊗(p::P{T,X}, q::P{S,Y}) is specic to basis/storage


basistype(p::AbstractUnivariatePolynomial{B,T,X}) where {B,T,X} = B
Base.eltype(p::AbstractUnivariatePolynomial{B,T,X}) where {B,T,X} = T
indeterminate(p::P) where {P <: AbstractUnivariatePolynomial} = indeterminate(P)
_indeterminate(::Type{P}) where {P <: AbstractUnivariatePolynomial} = nothing
_indeterminate(::Type{P}) where {B,T, X, P <: AbstractUnivariatePolynomial{B,T,X}} = X
indeterminate(::Type{P}) where {P <: AbstractUnivariatePolynomial} = something(_indeterminate(P), :x)

constructorof(::Type{<:AbstractUnivariatePolynomial}) = XXX()
⟒(P::Type{<:AbstractUnivariatePolynomial})  = constructorof(P) # returns the Storage{Basis} partially constructed type
⟒(p::P) where {P <: AbstractUnivariatePolynomial} = ⟒(P)

## Julia generics treating coefficients as an abstract vector
# pairs should iterate i => cᵢ where basis(P,i) is the basis vector
# * may not be in increasing or decreasing i
# * for standardbasis i -> xⁱ
# * possibly skipping when iszero(cᵢ)
Base.keys(p::AbstractUnivariatePolynomial) = Base.Generator(first, pairs(p))
Base.values(p::AbstractUnivariatePolynomial) = Base.Generator(last, pairs(p))
Base.firstindex(p::AbstractUnivariatePolynomial{B, T, X}) where {B,T,X} = XXX()
Base.lastindex(p::AbstractUnivariatePolynomial{B, T, X}) where {B,T,X} = XXX()
Base.iterate(p::AbstractUnivariatePolynomial, args...) = Base.iterate(p.coeffs, args...)
Base.pairs(p::AbstractUnivariatePolynomial) = XXX()

#Base.eltype(::Type{<:AbstractUnivariatePolynomial}) = Float64
Base.eltype(::Type{<:AbstractUnivariatePolynomial{B,T}}) where {B,T} = T

Base.size(p::AbstractUnivariatePolynomial) = (length(p),)
Base.size(p::AbstractUnivariatePolynomial, i::Integer) =  i <= 1 ? size(p)[i] : 1

Base.copy(p::AbstractUnivariatePolynomial) = XXX()

# dense collection
Base.iterate(p::AbstractUnivariatePolynomial, state = firstindex(p)) = _iterate(p, state) # _iterate in common.jl


#hasnan(p::AbstractUnivariatePolynomial) = any(hasnan, p)

# norm(q1 - q2)
function normΔ(q1::AbstractUnivariatePolynomial, q2::AbstractUnivariatePolynomial, p::Real = 2)
    iszero(q1) && return norm(q2, p)
    iszero(q2) && return norm(q1, p)
    r = zero(q1[end] + q2[end])
    tot = zero(r)
    for i ∈ union(keys(q1), keys(q2))
       @inbounds tot += (q1[i] - q2[i])^p
    end
    return tot^(1/p)
end


# map Polynomial terms -> vector terms
degree(p::AbstractUnivariatePolynomial) = iszero(p) ? -1 : lastindex(p)
# order(p::AbstractUnivariatePolynomial) = firstindex(p) XXX conflicts with DataFrames.order

# this helps, along with _set, make some storage-generic methods
_zeros(p::P, z, N) where {P <: AbstractUnivariatePolynomial} = _zeros(P, z, N)
_set(c::Vector, i, val)  = (c[i] = val; c)
_set(c::AbstractDict, i, val)  = (c[i] = val; c)
function _set(c::Tuple, i, val)
    @set! c[i] = val
    c
end

#check_same_variable(p::AbstractUnivariatePolynomial, q::AbstractUnivariatePolynomial) = indeterminate(p) == indeterminate(q)

# The zero polynomial. Typically has no coefficients
#Base.zero(p::P,args...) where {P <: AbstractUnivariatePolynomial} = zero(P,args...)
Base.zero(::Type{P}) where {B,P <: AbstractUnivariatePolynomial{B}} = zero(⟒(P){eltype(P),indeterminate(P)})
Base.zero(::Type{P},var::SymbolLike) where {B,P <: AbstractUnivariatePolynomial{B}} = zero(⟒(P){eltype(P),Symbol(var)})

# the polynomial 1
#Base.one(p::P,args...) where {P <: AbstractUnivariatePolynomial} = one(P,args...)
Base.one(::Type{P}) where {B, P <: AbstractUnivariatePolynomial{B}} = one(⟒(P){eltype(P),indeterminate(P)})
Base.one(::Type{P}, var::SymbolLike) where {B, P <: AbstractUnivariatePolynomial{B}} = one(⟒(P){eltype(P),Symbol(var)})

# the variable x
#variable(p::P) where {P <: AbstractUnivariatePolynomial} = variable(P)
variable(::Type{P}) where {B,P <: AbstractUnivariatePolynomial{B}} = variable(⟒(P){eltype(P),indeterminate(P)})
variable(::Type{P}, var::SymbolLike) where {B,P<:AbstractUnivariatePolynomial{B}} = variable(⟒(P){eltype(P),Var(var)})

# i -> basis polynomial
basis(p::P, i::Int) where {P <: AbstractUnivariatePolynomial} = basis(P, i)
basis(::Type{P}, i::Int) where {B,P <: AbstractUnivariatePolynomial{B}} = basis(⟒(P){eltype(P),indeterminate(P)}, i)

_convert(p::P, as) where {B,T,X,P <: AbstractUnivariatePolynomial{B,T,X}} = ⟒(P){T,X}(as)
copy_with_eltype(::Type{T}, ::Val{X}, p::P) where {B,T, X, S, Y, P <:AbstractUnivariatePolynomial{B,S,Y}} =
    ⟒(P){T, Symbol(X)}(p.coeffs)

# return dense coefficients (vector or tuple)
coeffs(p::AbstractUnivariatePolynomial) = [p[i] for i ∈ firstindex(p):lastindex(p)]

# function isconstant(p::AbstractUnivariatePolynomial)
#     p₀ = trim_trailing_zeros(p)
#     return (firstindex(p₀) == lastindex(p₀) == 0)
# end

# chop chops right side of p
# use trunc for left and right
# can pass tolerances
Base.chop(p::AbstractUnivariatePolynomial; kwargs...) = chop!(copy(p))
chop!(p::AbstractUnivariatePolynomial; kwargs...) = XXX()

## --- constant term ---

# arithmetic dispatch
struct ConstantTerm{T}
    x::T
end
function ConstantTerm(p::AbstractUnivariatePolynomial)
    isconstant(p) || throw(ArgumentError("Non-constant polynomial"))
    convert(ConstantTerm, p)
end
Base.getindex(b::ConstantTerm) = b.x
Base.show(io::IO, c::ConstantTerm) = print(io, c.x)
Base.iszero(b::ConstantTerm) = iszero(b.x)
isconstant(::ConstantTerm) = true
Base.convert(::Type{ConstantTerm}, p::AbstractUnivariatePolynomial) = ConstantTerm(constantterm(p))
Base.convert(::Type{ConstantTerm{T}}, p::AbstractUnivariatePolynomial) where {T} = ConstantTerm(T(constantterm(p)))


## ---

#= Comparisons =#
# need to promote Number -> Poly
# Base.isapprox(p1::AbstractUnivariatePolynomial, p2::Number; kwargs...) = isapprox(promote(p1, p2)...; kwargs...)
# Base.isapprox(p1::Number, p2::AbstractUnivariatePolynomial; kwargs...) = isapprox(promote(p1, p2)...; kwargs...)
function Base.isapprox(p1::AbstractUnivariatePolynomial, p2::AbstractUnivariatePolynomial; kwargs...)
    isapprox(promote(p1, p2)...; kwargs...)
end

function Base.isapprox(p1::AbstractUnivariatePolynomial{B}, p2::AbstractUnivariatePolynomial{B}; kwargs...) where {B}
    if isconstant(p1)
        isconstant(p2) && return constantterm(p1) == constantterm(p2)
        return false
    elseif isconstant(p2)
        return false
    end
    assert_same_variable(p1, p2) || return false
    isapprox(promote(p1, p2)...; kwargs...)
end

function Base.isapprox(p1::AbstractUnivariatePolynomial{B,T,X},
                       p2::AbstractUnivariatePolynomial{B,T,X};
                       rtol::Real = (Base.rtoldefault(T,T,0)),
                       atol::Real = 0,) where {B,T,X}
    (hasnan(p1) || hasnan(p2)) && return false  # NaN poisons comparisons
    # copy over from abstractarray.jl
    Δ  = normΔ(p1,p2)
    if isfinite(Δ)
        return Δ <= max(atol, rtol * max(norm(p1), norm(p2)))
    else
        for i in 0:max(degree(p1), degree(p2))
            isapprox(p1[i], p2[i]; rtol=rtol, atol=atol) || return false
        end
        return true
    end
end

function Base.isapprox(p1::AbstractUnivariatePolynomial{B,T,X}, p2::Number; kwargs...) where {B,T,X}
    q = p2 * one(⟒(p1){T,X})
    isapprox(p1, q; kwargs...)
end
Base.isapprox(p1::Number, p2::AbstractUnivariatePolynomial; kwargs...) = isapprox(p2, p1; kwargs...)

Base.isequal(p1::P, p2::P) where {P <: AbstractUnivariatePolynomial} = hash(p1) == hash(p2)
function Base.:(==)(p1::P, p2::P) where {P <: AbstractUnivariatePolynomial}
    iszero(p1) && iszero(p2) && return true
    lastindex(p1) == lastindex(p2) || return false
    # coeffs(p1) == coeffs(p2), but non-allocating
    for i ∈ union(keys(p1), keys(p2))
        p1[i] == p2[i] || return false
#    for ((i,pᵢ), (j, pⱼ)) ∈ zip(pairs(p1), pairs(p2))
#        i == j && pᵢ == pⱼ || return false
    end
    return true
end
function Base.:(==)(p1::AbstractUnivariatePolynomial, p2::AbstractUnivariatePolynomial)
    if isconstant(p1)
        isconstant(p2) && return constantterm(p1) == constantterm(p2)
        return false
    elseif isconstant(p2)
        return false # p1 is not constant
    end
    check_same_variable(p1, p2) || return false
    ==(promote(p1,p2)...)
end
Base.:(==)(p::AbstractUnivariatePolynomial, n::Number) = degree(p) <= 0 && constantterm(p) == n
Base.:(==)(n::Number, p::AbstractUnivariatePolynomial) = p == n

## --- arithmetic operations ---
## implement
## * unary - : here using scalar_mutl
## * scalar_add : with basis
## * scalar_mult : with storage type
## * scalar division: here using scalar_mult
## * polynomial addition: with storage type
## * polynomial multiplication: resolstorage type + basis
##
Base.:-(p::AbstractUnivariatePolynomial) = scalar_mult(-1, p)

Base.:+(c::Scalar, p::AbstractUnivariatePolynomial) = scalar_add(p, c)
Base.:+(p::AbstractUnivariatePolynomial, c::Scalar) = scalar_add(p, c)
Base.:+(c::ConstantTerm, p::AbstractUnivariatePolynomial) = scalar_add(p, c[])
Base.:+(p::AbstractUnivariatePolynomial, c::ConstantTerm) = scalar_add(p, c[])
scalar_add(p::AbstractUnivariatePolynomial, c) = scalar_add(c,p) # scalar addition is commutative

Base.:+(p::AbstractUnivariatePolynomial) = p
Base.:+(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, X}) where {B,T,S,X} =
            +(promote(p,q)...)
Base.:+(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, Y}) where {B,T,S,X,Y} =
            _mixed_symbol_op(+, p, q)

Base.:-(c::Scalar, p::AbstractUnivariatePolynomial) = c + (-p)
Base.:-(p::AbstractUnivariatePolynomial, c::Scalar) = p + (-c)
Base.:-(c::ConstantTerm, p::AbstractUnivariatePolynomial) = (-c[]) + p
Base.:-(p::AbstractUnivariatePolynomial, c::ConstantTerm) = p - c[]

Base.:-(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, X}) where {B,T,S,X} = p + (-q)
Base.:-(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, Y}) where {B,T,S,X,Y} =
            _mixed_symbol_op(-, p, q)

Base.:*(c::Scalar, p::ConstantTerm) = ConstantTerm(c*p[])
Base.:*(p::ConstantTerm, c::Scalar) = ConstantTerm(p[] * c)

Base.:*(c::Scalar, p::AbstractUnivariatePolynomial) = scalar_mult(c, p)
Base.:*(c::ConstantTerm, p::AbstractUnivariatePolynomial) = scalar_mult(c[], p)
Base.:*(p::AbstractUnivariatePolynomial, c::Scalar) = scalar_mult(p, c)
Base.:*(p::AbstractUnivariatePolynomial, c::ConstantTerm) = scalar_mult(p, c[])

Base.:*(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, X}) where {B,T,S,X} = *(promote(p,q)...)
Base.:*(p::P, q::P) where {B,T,X,P <: AbstractUnivariatePolynomial{B,T,X}} =
            ⊗(p, q)
Base.:*(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, Y}) where {B,T,S,X,Y}  =
            _mixed_symbol_op(*, p, q)

Base.:/(p::AbstractUnivariatePolynomial, c::Scalar) = scalar_mult(p, one(eltype(p))/c)
Base.:/(p::AbstractUnivariatePolynomial, c::ConstantTerm) = scalar_mult(p, one(eltype(p))/c)

Base.:^(p::AbstractUnivariatePolynomial, n::Integer) = Base.power_by_squaring(p, n)

# treat constant polynomials as ConstantTerm when symbols mixed
function _mixed_symbol_op(op,
                          p::AbstractUnivariatePolynomial{B, T, X},
                          q::AbstractUnivariatePolynomial{B, S, Y}) where {B,T,S,X,Y}
    X == Y && throw(ArgumentError("dispatch should catch this case"))
    if isconstant(p)
        return  op(convert(ConstantTerm, p), q)
    elseif isconstant(q)
        return  op(p, convert(ConstantTerm, q))
    end
    assert_same_variable(X,Y)
end


# only need to define differentiate(p::PolyType)
function derivative(p::AbstractUnivariatePolynomial, n::Int=1)
    n < 0 && throw(ArgumentError("n must be non-negative"))
    iszero(n) && return p
    p′ = differentiate(p)
    for i ∈ 2:n
        p′ = differentiate(p′)
    end
    p′
end
const differentiate = derivative

# only need to define integrate(p::PolyType)
function integrate(p::AbstractUnivariatePolynomial, c)
    scalar_add(integrate(p), c)
end





macro poly_register(name)
    poly = esc(name)
    quote
        #Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote(p::P, q::Q) where {B, X, T, P <:$poly{B,T,X}, Q <: $poly{B,T,X}} = p,q
        Base.promote_rule(::Type{<:$poly{B,T,X}}, ::Type{<:$poly{B,S,X}}) where {B,T,S,X} =  $poly{B,promote_type(T, S),X}
        Base.promote_rule(::Type{<:$poly{B,T,X}}, ::Type{S}) where {B,T,S<:Scalar,X} =
            $poly{B,promote_type(T, S), X}

        # constants
        $poly{B,T,X}(n::S) where {B, T, X, S<:Scalar} =
            T(n) *  one($poly{B, T, X})
        $poly{B, T}(n::S, var::SymbolLike = Var(:x)) where {B, T, S<:Scalar} =
            T(n) *  one($poly{B, T, Symbol(var)})
        $poly{B}(n::S, var::SymbolLike = Var(:x))  where {B, S  <: Scalar} = n * one($poly{B, S, Symbol(var)})

        $poly{B,T}(var::SymbolLike=Var(:x)) where {B,T} = variable($poly{B, T, Symbol(var)})
        $poly{B}(var::SymbolLike=Var(:x)) where {B} = variable($poly{B}, Symbol(var))

        (p::$poly)(x) = evalpoly(x, p)
    end
end
