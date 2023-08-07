"""
    AbstractUnivariatePolynomial{B,T,X} <: AbstractPolynomial{T,X}
    AbstractDenseUnivariatePolynomial{B,T,X} <: AbstractUnivariatePolynomial{B,T,X}
    AbstractLaurentUnivariatePolynomial{B,T,X} <: AbstractUnivariatePolynomial{B,T,X}

Abstract container types for polynomials with an explicit basis, `B`.
`AbstractDenseUnivariatePolynomial` is for `0`-based polynomials;
`AbstractLaurentUnivariatePolynomial` is for polynomials with possibly negative powers of the indeterminate.

"""
abstract type AbstractUnivariatePolynomial{B, T, X} <: AbstractPolynomial{T,X} end

# for 0-based polys
abstract type AbstractDenseUnivariatePolynomial{B, T, X} <: AbstractUnivariatePolynomial{B,T,X} end

# for negative integer powers
abstract type AbstractLaurentUnivariatePolynomial{B, T, X} <: AbstractUnivariatePolynomial{B,T,X} end

"""
    AbstractBasis

Abstract type for specifying a polynomial basis.
"""
abstract type AbstractBasis end
export AbstractUnivariatePolynomial

XXX() = throw(ArgumentError("No default method defined"))

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

# type of basis
basistype(p::AbstractUnivariatePolynomial{B,T,X}) where {B,T,X} = B
basistype(::Type{<:AbstractUnivariatePolynomial{B}}) where {B} = B # some default

# eltype of polynomial
Base.eltype(p::AbstractUnivariatePolynomial{B,T,X}) where {B,T,X} = T

# indeterminate of polynomial
indeterminate(p::P) where {P <: AbstractUnivariatePolynomial} = indeterminate(P)
_indeterminate(::Type{P}) where {P <: AbstractUnivariatePolynomial} = nothing
_indeterminate(::Type{P}) where {B,T, X, P <: AbstractUnivariatePolynomial{B,T,X}} = X
indeterminate(::Type{P}) where {P <: AbstractUnivariatePolynomial} = something(_indeterminate(P), :x)


constructorof(::Type{<:AbstractUnivariatePolynomial}) = XXX() # defined in container-types
⟒(P::Type{<:AbstractUnivariatePolynomial})  = constructorof(P)
⟒(p::P) where {P <: AbstractUnivariatePolynomial} = ⟒(P)

## Julia generics treating coefficients as an abstract vector
# pairs should iterate i => cᵢ where basis(P,i) is the basis vector
# * may not be in increasing or decreasing i
# * for standardbasis i -> xⁱ
# * possibly skipping when iszero(cᵢ)
#Base.firstindex(p::AbstractUnivariatePolynomial{B, T, X}) where {B,T,X} = XXX()
#Base.lastindex(p::AbstractUnivariatePolynomial{B, T, X}) where {B,T,X} = XXX()
#Base.iterate(p::AbstractUnivariatePolynomial, args...) = Base.iterate(values(p), args...)
Base.iterate(p::AbstractUnivariatePolynomial, state = firstindex(p)) = _iterate(p, state) # _iterate in common.jl

Base.keys(p::AbstractUnivariatePolynomial)   = eachindex(p)
Base.values(p::AbstractUnivariatePolynomial) = values(p.coeffs) # Dict based containers must specialize
Base.pairs(p::AbstractUnivariatePolynomial)  = Base.Generator(=>, keys(p), values(p))



#Base.eltype(::Type{<:AbstractUnivariatePolynomial}) = Float64 # common.jl
Base.eltype(::Type{<:AbstractUnivariatePolynomial{B,T}}) where {B,T} = T

Base.size(p::AbstractUnivariatePolynomial) = (length(p),)
Base.size(p::AbstractUnivariatePolynomial, i::Integer) =  i <= 1 ? size(p)[i] : 1

#Base.copy(p::AbstractUnivariatePolynomial) = XXX()

# map Polynomial terms -> vector terms
# Default degree **assumes** basis element Tᵢ has degree i.
degree(p::AbstractUnivariatePolynomial) = iszero(p) ? -1 : lastindex(p) # XXX() is likely a safer choice

# this helps, along with _set, make some storage-generic methods
_zeros(p::P, z, N) where {P <: AbstractUnivariatePolynomial} = _zeros(P, z, N)
_zeros(::Type{<:AbstractDenseUnivariatePolynomial}, z::T, N) where {T} = fill(z,(N,)) # default
_set(c::Vector, i, val)  = (c[i] = val; c)
_set(c::AbstractDict, i, val)  = (c[i] = val; c)
function _set(c::Tuple, i, val)
    @set! c[i] = val
    c
end

# for indexing for Tᵢ i + offset should be array index
# this default is for 1-based backends (vector, tuple)
offset(p::AbstractDenseUnivariatePolynomial) = 1

# i -> basis polynomial. Uses `order` argument, which may save some space
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
#chop!(p::AbstractUnivariatePolynomial; kwargs...) = XXX()
chop!(p::AbstractDenseUnivariatePolynomial; kwargs...) = (chop!(p.coeffs); p) # default for mutable vector backed; tuple backed need other



## ---

#= Comparisons =#
# iterable over keys of both
function keys_union(p::AbstractUnivariatePolynomial, q::AbstractUnivariatePolynomial)
    minimum(firstindex, (p,q)):maximum(lastindex, (p,q))
end


# norm(q1 - q2) only non-allocating (save for unique)
function normΔ(q1::AbstractUnivariatePolynomial, q2::AbstractUnivariatePolynomial)
    iszero(q1) && return norm(q2, 2)
    iszero(q2) && return norm(q1, 2)
    tot = abs(zero(q1[end] + q2[end]))
    for i ∈ keys_union(q1, q2)
       @inbounds tot += abs2(q1[i] - q2[i])
    end
    return sqrt(tot)
end

# need to promote Number -> Poly
function Base.isapprox(p1::AbstractUnivariatePolynomial, p2::AbstractUnivariatePolynomial; kwargs...)
    isapprox(promote(p1, p2)...; kwargs...)
end

# compare X != Y
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
                       p2::AbstractUnivariatePolynomial{B,S,X};
                       rtol = nothing,
                       atol = nothing) where {B,T,S, X}

    (hasnan(p1) || hasnan(p2)) && return false  # NaN poisons comparisons
    R = float(real(promote_type(T,S)))
    rtol = something(rtol, Base.rtoldefault(R,R,0))
    atol = something(atol, 0)


    # copy over from abstractarray.jl
    Δ  = normΔ(p1,p2)
    if isfinite(Δ)
        return Δ <= max(atol, rtol * max(norm(p1), norm(p2)))
    else
        for i in keys_union(p1, p2)
            isapprox(p1[i], p2[i]; rtol=rtol, atol=atol) || return false
        end
        return true
    end
end


## --- arithmetic operations ---
## all in common.jl

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
## better parallel with integrate, but derivative has been used in this package for too long to change
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
# Base.one(::Type{P}) where {B,T,X,P<:AbstractUnivariatePolynomial{B,T,X}} = XXX()
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
        Polynomials.constructorof(::Type{<:$poly{B}}) where {B} = $poly{B}
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
