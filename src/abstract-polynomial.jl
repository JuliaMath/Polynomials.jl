# XXX Maybe merge into with common.jl, abstract.jl
"""
    AbstractUnivariatePolynomial{B,T,X} <: AbstractPolynomial{T,X}

Abstract type for polynomials with an explicit basis, `B`.
"""
abstract type AbstractUnivariatePolynomial{B, T, X} <: AbstractPolynomial{T,X} end

# for 0-based polys
abstract type AbstractDenseUnivariatePolynomial{B, T, X} <: AbstractUnivariatePolynomial{B,T,X} end

# for negative integer powers
abstract type AbstractLaurentUnivariatePolynomial{B, T, X} <: AbstractUnivariatePolynomial{B,T,X} end

abstract type AbstractBasis end
export AbstractUnivariatePolynomial

XXX() = throw(ArgumentError("Method not defined"))

## --------------------------------------------------

function showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {B, T, P<:AbstractUnivariatePolynomial{B,T}}
    if _iszero(pj) return false end

    pj = printsign(io, pj, first, mimetype)
    if hasone(T)
        if !(_isone(pj) && !(showone(T) || j == 0))
            printcoefficient(io, pj, j, mimetype)
        end
    else
        printcoefficient(io, pj, j, mimetype)
    end

    printproductsign(io, pj, j, mimetype)
    #printbasis(io, P, j, mimetype)
    printbasis(io, ⟒(P){T,Symbol(var)}, j, mimetype)
    return true
end

# overload printbasis to see a subscript, as Tᵢ... or if there are parameters in basis
function printbasis(io::IO, ::Type{P}, j::Int, m::MIME) where {P <: AbstractUnivariatePolynomial}
    print(io, basis_symbol(P))
    print(io, subscript_text(j, m))
end

basis_symbol(::Type{AbstractUnivariatePolynomial{B,T,X}}) where {B,T,X} = "Χ($(X))"


## idea is vector space stuff (scalar_add, scalar_mult, vector +/-, ^) goes here
## connection (convert, transform) is specific to a basis (storage)
## ⊗(p::P{T,X}, q::P{S,Y}) is specic to basis/storage


basistype(p::AbstractUnivariatePolynomial{B,T,X}) where {B,T,X} = B
basistype(::Type{<:AbstractUnivariatePolynomial{B}}) where {B} = B # some default
Base.eltype(p::AbstractUnivariatePolynomial{B,T,X}) where {B,T,X} = T
indeterminate(p::P) where {P <: AbstractUnivariatePolynomial} = indeterminate(P)
_indeterminate(::Type{P}) where {P <: AbstractUnivariatePolynomial} = nothing
_indeterminate(::Type{P}) where {B,T, X, P <: AbstractUnivariatePolynomial{B,T,X}} = X
indeterminate(::Type{P}) where {P <: AbstractUnivariatePolynomial} = something(_indeterminate(P), :x)


constructorof(::Type{<:AbstractUnivariatePolynomial}) = XXX()
⟒(P::Type{<:AbstractUnivariatePolynomial})  = constructorof(P)
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
#Base.iterate(p::AbstractUnivariatePolynomial, args...) = Base.iterate(values(p), args...)
Base.iterate(p::AbstractUnivariatePolynomial, state = firstindex(p)) = _iterate(p, state) # _iterate in common.jl
Base.pairs(p::AbstractUnivariatePolynomial) = XXX()

#Base.eltype(::Type{<:AbstractUnivariatePolynomial}) = Float64
Base.eltype(::Type{<:AbstractUnivariatePolynomial{B,T}}) where {B,T} = T

Base.size(p::AbstractUnivariatePolynomial) = (length(p),)
Base.size(p::AbstractUnivariatePolynomial, i::Integer) =  i <= 1 ? size(p)[i] : 1

Base.copy(p::AbstractUnivariatePolynomial) = XXX()

#hasnan(p::AbstractUnivariatePolynomial) = any(hasnan, p)

# map Polynomial terms -> vector terms
# Default degree **assumes** basis element Tᵢ has degree i.
degree(p::AbstractUnivariatePolynomial) = iszero(p) ? -1 : lastindex(p) # XXX() is likely a safer choice

# this helps, along with _set, make some storage-generic methods
_zeros(p::P, z, N) where {P <: AbstractUnivariatePolynomial} = _zeros(P, z, N)
_set(c::Vector, i, val)  = (c[i] = val; c)
_set(c::AbstractDict, i, val)  = (c[i] = val; c)
function _set(c::Tuple, i, val)
    @set! c[i] = val
    c
end

#check_same_variable(p::AbstractUnivariatePolynomial, q::AbstractUnivariatePolynomial) = indeterminate(p) == indeterminate(q)

# The zero polynomial. Typically has no coefficients. in common.jl
#Base.zero(p::P,args...) where {P <: AbstractUnivariatePolynomial} = zero(P,args...)
#Base.zero(::Type{P}) where {B,P <: AbstractUnivariatePolynomial{B}} = zero(⟒(P){eltype(P),indeterminate(P)})
#Base.zero(::Type{P},var::SymbolLike) where {B,P <: AbstractUnivariatePolynomial{B}} = zero(⟒(P){eltype(P),Symbol(var)})


# the polynomial 1
# one(P) is basis dependent and must be implemented in one(::Type{<:P})
Base.one(p::P,args...) where {P <: AbstractUnivariatePolynomial} = one(P,args...)
Base.one(::Type{P}) where {B, P <: AbstractUnivariatePolynomial{B}} = one(⟒(P){eltype(P),indeterminate(P)})
Base.one(::Type{P}, var::SymbolLike) where {B, P <: AbstractUnivariatePolynomial{B}} = one(⟒(P){eltype(P),Symbol(var)})

# the variable x is basid dependent and must be implmented in variable(::Type{<:P})
variable(p::P) where {P <: AbstractUnivariatePolynomial} = variable(P)
variable(::Type{P}) where {B,P <: AbstractUnivariatePolynomial{B}} = variable(⟒(P){eltype(P),indeterminate(P)})
variable(::Type{P}, var::SymbolLike) where {B,P<:AbstractUnivariatePolynomial{B}} = variable(⟒(P){eltype(P),Var(var)})

# i -> basis polynomial
basis(p::P, i::Int) where {P <: AbstractUnivariatePolynomial} = basis(P, i)
basis(::Type{P}, i::Int) where {B,P <: AbstractUnivariatePolynomial{B}} = basis(⟒(P){eltype(P),indeterminate(P)}, i)

copy_with_eltype(::Type{T}, ::Val{X}, p::P) where {B,T, X, S, Y, P <:AbstractUnivariatePolynomial{B,S,Y}} =
    ⟒(P){T, Symbol(X)}(p.coeffs)


coeffs(p::P) where {P <: AbstractLaurentUnivariatePolynomial} = p.coeffs
function coeffs(p::P) where {P <: AbstractDenseUnivariatePolynomial}
    firstindex(p) == 0 && return p.coeffs
    firstindex(p) > 0 && return [p[i] for i ∈ 0:lastindex(p)]
    throw(ArgumentError("Polynomial type does not support negative degree terms"))
end

gtτ(x, τ) = abs(x) > τ
# return index or nothing of last non "zdero"
function chop_right_index(x; rtol=nothing, atol=nothing)
    isempty(x) && return nothing
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, norm(x,2) * δ)
    i = findlast(Base.Fix2(gtτ, τ), x)
    i
end

# chop chops right side of p
# use trunc for left and right
# can pass tolerances
Base.chop(p::AbstractUnivariatePolynomial; kwargs...) = chop!(copy(p))
chop!(p::AbstractUnivariatePolynomial; kwargs...) = XXX()



## ---

#= Comparisons =#

# norm(q1 - q2) only non-allocating
function normΔ(q1::AbstractUnivariatePolynomial, q2::AbstractUnivariatePolynomial)
    iszero(q1) && return norm(q2, 2)
    iszero(q2) && return norm(q1, 2)
    tot = abs(zero(q1[end] + q2[end]))
    for i ∈ unique(Base.Iterators.flatten((keys(q1), keys(q2))))
       @inbounds tot += abs2(q1[i] - q2[i])
    end
    return sqrt(tot)
end

# need to promote Number -> Poly
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


## --- arithmetic operations ---
## implement
## * unary - : here using scalar_mutl
## * scalar_add : with basis
## * scalar_mult : with storage type
## * scalar division: here using scalar_mult
## * polynomial addition: with storage type
## * polynomial multiplication: resolstorage type + basis
##
Base.:-(p::AbstractUnivariatePolynomial) = map(-, p) #scalar_mult(-1, p)

Base.:+(c::Scalar, p::AbstractUnivariatePolynomial) = scalar_add(c, p)
Base.:+(p::AbstractUnivariatePolynomial, c::Scalar) = scalar_add(c, p)
scalar_add(p::AbstractUnivariatePolynomial, c) = scalar_add(c, p) # scalar addition is commutative

Base.:+(p::AbstractUnivariatePolynomial) = p
Base.:+(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, X}) where {B,T,S,X} =
            +(promote(p,q)...)
Base.:+(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, Y}) where {B,T,S,X,Y} =
            _mixed_symbol_op(+, p, q)

Base.:-(c::Scalar, p::AbstractUnivariatePolynomial) = c + (-p)
Base.:-(p::AbstractUnivariatePolynomial, c::Scalar) = p + (-c)

Base.:-(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, X}) where {B,T,S,X} = p + (-q)
Base.:-(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, Y}) where {B,T,S,X,Y} =
            _mixed_symbol_op(-, p, q)

Base.:*(c::Scalar, p::AbstractUnivariatePolynomial) = scalar_mult(c, p)
Base.:*(p::AbstractUnivariatePolynomial, c::Scalar) = scalar_mult(p, c)

Base.:*(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, X}) where {B,T,S,X} = *(promote(p,q)...)
Base.:*(p::P, q::P) where {B,T,X,P <: AbstractUnivariatePolynomial{B,T,X}} =
            ⊗(p, q)
Base.:*(p::AbstractUnivariatePolynomial{B, T, X},
        q::AbstractUnivariatePolynomial{B, S, Y}) where {B,T,S,X,Y}  =
            _mixed_symbol_op(*, p, q)

Base.:/(p::AbstractUnivariatePolynomial, c::Scalar) = scalar_div(p, c)

Base.:^(p::AbstractUnivariatePolynomial, n::Integer) = Base.power_by_squaring(p, n)

scalar_mult(p::AbstractUnivariatePolynomial{B,T,X}, c) where {B,T,X} = map(Base.Fix2(*,c), p)
scalar_mult(c, p::AbstractUnivariatePolynomial{B,T,X}) where {B,T,X} = map(Base.Fix1(*,c), p)
scalar_div(p::AbstractUnivariatePolynomial{B,T,X}, c) where {B,T,X}  = map(Base.Fix2(*,one(T)/c), p) # much better than Fix2(/,c)

# treat constant polynomials as constants when symbols mixed
scalar_op(::typeof(*)) = scalar_mult
scalar_op(::typeof(+)) = scalar_add
scalar_op(::typeof(/)) = scalar_div
function _mixed_symbol_op(op,
                          p::P,
                          q::Q) where {B,T,S,X,Y,
                                       P<:AbstractUnivariatePolynomial{B, T, X},
                                       Q<:AbstractUnivariatePolynomial{B, S, Y}}
    X == Y && throw(ArgumentError("dispatch should catch this case"))
    isconstant(p) && return scalar_op(op)(constantterm(p), q)
    isconstant(q) && return scalar_op(op)(p, constantterm(q))
    assert_same_variable(X,Y)
end


# only need to define derivative(p::PolyType)
function derivative(p::AbstractUnivariatePolynomial, n::Int=1)
    n < 0 && throw(ArgumentError("n must be non-negative"))
    iszero(n) && return p
    p′ = differentiate(p)
    for i ∈ 2:n
        p′ = differentiate(p′)
    end
    p′
end
## better parallel with integrate, but derivative has been used here.
const differentiate = derivative

function fit(::Type{P},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg,
             cs::Dict;
             kwargs...) where {T, P<:AbstractUnivariatePolynomial}
    convert(P, fit(Polynomial, x, y, deg, cs; kwargs...))
end


## Interface
## These must be implemented for a storage type / basis
# minimumexponent(::Type{P}) where {B,P<:AbstractUnivariatePolynomial{B}} = XXX()
# Base.one(::Type{P}) where {B,P<:AbstractUnivariatePolynomial{B}} = XXX()
# variable(::Type{P}) where {B,P<:AbstractUnivariatePolynomial{B}} = XXX()
# evalpoly(x, p::P) where {B,P<:AbstractUnivariatePolynomial{B}} = XXX()
# scalar_add(c, p::P) where {B,P<:AbstractUnivariatePolynomial{B}} = XXX()
# ⊗(p::P, q::Q) where {B,P<:AbstractUnivariatePolynomial{B},Q<:AbstractUnivariatePolynomial{B}} = XXX()

# these *may* be implemented for a basis type
# * basis_symbol/printbasis
# * one, variable, constantterm, domain, mapdomain
# * derivative
# * integrate
# * divrem
# * vander
# * companion

# promote, promote_rule, vector specification, untyped specification, handle constants,  conversion of Q(p)
# poly composition, calling a polynomial
macro poly_register(name)
    poly = esc(name)
    quote
        #Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote(p::P, q::Q) where {B, X, T, P <:$poly{B,T,X}, Q <: $poly{B,T,X}} = p,q
        Base.promote_rule(::Type{<:$poly{B,T,X}}, ::Type{<:$poly{B,S,X}}) where {B,T,S,X} =  $poly{B,promote_type(T, S),X}
        Base.promote_rule(::Type{<:$poly{B,T,X}}, ::Type{S}) where {B,T,S<:Scalar,X} =
            $poly{B,promote_type(T, S), X}

        # vector
        $poly{B,T}(xs::AbstractVector{S}, order::Int, var::SymbolLike=Var(:x)) where {B,T,S} = $poly{B,T,Symbol(var)}(xs,order)
        $poly{B,T}(xs::AbstractVector{S}, var::SymbolLike) where {B,T,S} = $poly{B,T,Symbol(var)}(xs,0)
        $poly{B}(xs::AbstractVector{T}, order::Int, var::SymbolLike=Var(:x)) where {B,T} = $poly{B,T,Symbol(var)}(xs,order)
        $poly{B}(xs::AbstractVector{T}, var::SymbolLike) where {B,T} = $poly{B,T,Symbol(var)}(xs,0)

        # untyped
        $poly{B,T,X}(xs, order::Int=0) where {B,T,X} = $poly{B,T,X}(collect(T,xs), order)
        $poly{B,T}(xs, order::Int=0, var::SymbolLike=Var(:x)) where {B,T} = $poly{B,T,Var(var)}(collect(T,xs), order)
        $poly{B,T}(xs, var::SymbolLike) where {B,T} = $poly{B,T}(xs, 0, var)
        function $poly{B}(xs, order::Int=0, var::SymbolLike=Var(:x)) where {B}
            cs = collect(promote(xs...))
            T = eltype(cs)
            $poly{B,T,Symbol(var)}(cs, order)
        end
        $poly{B}(xs, var::SymbolLike) where {B} = $poly{B}(xs, 0, var)

        # constants
        $poly{B,T,X}(n::S) where {B, T, X, S<:Scalar} =
            T(n) *  one($poly{B, T, X})
        $poly{B, T}(n::S, var::SymbolLike = Var(:x)) where {B, T, S<:Scalar} =
            T(n) *  one($poly{B, T, Symbol(var)})
        $poly{B}(n::S, var::SymbolLike = Var(:x))  where {B, S  <: Scalar} = n * one($poly{B, S, Symbol(var)})

        $poly{B,T}(var::SymbolLike=Var(:x)) where {B,T} = variable($poly{B, T, Symbol(var)})
        $poly{B}(var::SymbolLike=Var(:x)) where {B} = variable($poly{B}, Symbol(var))

        # conversion via P(q)
        $poly{B,T,X}(c::AbstractPolynomial{S,Y}) where {B,T,X,S,Y} = convert($poly{B,T,X}, c)
        $poly{B,T}(c::AbstractPolynomial{S,Y}) where {B,T,S,Y} = convert($poly{B,T}, c)
        $poly{B}(c::AbstractPolynomial{S,Y}) where {B,S,Y} = convert($poly{B}, c)

        # poly composition and evaluation
        (p::$poly)(x::AbstractPolynomial) = polynomial_composition(p, x)
        (p::$poly)(x) = evalpoly(x, p)
    end
end
