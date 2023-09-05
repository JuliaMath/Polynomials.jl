export AbstractPolynomial

# *internal* means to pass variable symbol to constructor through 2nd position and keep type stability
struct Var{T} end
Var(x::Var) = x
Var(x::Symbol) = Var{x}()
Var(x::Type{Var{u}}) where {u} = x
Var(x::AbstractString) = Var(Symbol(x))
Var(x::Char) = Var(Symbol(x))

Base.Symbol(::Var{T}) where {T} = T

const SymbolLike = Union{AbstractString,Char,Symbol, Var{T} where T}


"""
    AbstractPolynomial{T,X}

An abstract type for various polynomials with an *implicit* basis.

A polynomial type holds an indeterminate `X`; coefficients of type `T`, stored in some container type; and an implicit basis, such as the standard basis.

# Properties
- `coeffs` - The coefficients of the polynomial

# The type `T`

`T` need not be `T <: Number`, at the moment it is not restricted

Some `T`s will not be successful

* scalar mult: `c::Number * p::Polynomial` should be defined
* scalar mult: `c::T * p::Polynomial{T}` An  ambiguity when `T <: AbstractPolynomial`
* scalar mult: `p::Polynomial{T} * c::T` need not commute

* add/sub: `p::Polynomial{T} + q::Polynomial{T}` should be defined
* sub: `p -p` sometimes needs `zero{T}` defined

* poly mult: `p::Polynomial{T} * q::Polynomial{T}` Needs "`T * T`" defined (e.g. `Base.promote_op(*, Vector{Int}, Vector{Int}))` needs to be something.)
* poly powers: `p::Polynomial{T}^2` needs "`T^2`" defined

* implicit promotion: `p::Polynomial{T} + c::Number`  needs `convert(T, c)` defined
* implicit promotion: `p::Polynomial{T} + c::T`  ambiguity if `T <: AbstractPolynomial`

* evaluation: `p::Polynomial{T}(s::Number)`
* evaluation `p::Polynomial{T}(c::T)`   needs `T*T` defined
* evaluation of a `0` polynomial requires `zero(T)` to be defined.

* derivatives: `derivative(p::Polynomial{T})`
* integrals: `integrate(p::Polynomial{T})`


"""
abstract type AbstractPolynomial{T,X} end


## -----
# We want  ⟒(P{α…,T,X}) = P{α…}; this default
# works for most cases
⟒(P::Type{<:AbstractPolynomial}) = constructorof(P)
⟒(p::P) where {P <: AbstractPolynomial} = ⟒(P)

"""
    Polynomials.@register(name)

Given a polynomial with `name`, creates some common convenience constructors and conversions to minimize code required for implementation of a new polynomial type.

# Example
```julia
struct MyPolynomial{T,X} <: AbstractPolynomial{T,X} end

Polynomials.@register MyPolynomial
```

# Implementations
This will implement simple self-conversions like `convert(::Type{MyPoly}, p::MyPoly) = p` and creates two promote rules. The first allows promotion between two types (e.g. `promote(Polynomial, ChebyshevT)`) and the second allows promotion between parameterized types (e.g. `promote(Polynomial{T}, Polynomial{S})`).

For constructors, it implements the shortcut for `MyPoly(...) = MyPoly{T}(...)`, singleton constructor `MyPoly(x::Number, ...)`,  conversion constructor `MyPoly{T}(n::S, ...)`, and `variable` alternative  `MyPoly(var=:x)`.
"""
macro register(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        function Base.convert(P::Type{<:$poly}, p::$poly{T,X}) where {T,X}
            isconstant(p) && return constructorof(P){eltype(P),indeterminate(P)}(constantterm(p))
            constructorof(P){eltype(P), indeterminate(P,p)}(coeffs(p))
        end
        Base.promote(p::P, q::Q) where {X, T, P <:$poly{T,X}, Q <: $poly{T,X}} = p,q
        Base.promote_rule(::Type{<:$poly{T,X}}, ::Type{<:$poly{S,X}}) where {T,S,X} =  $poly{promote_type(T, S),X}
        Base.promote_rule(::Type{P}, ::Type{S}) where {T,S<:Number,X,P<:$poly{T,X}} =
            $poly{promote_type(T, S),X}
#        Base.promote_rule(::Type{<:$poly{T,X}}, ::Type{S}) where {T,S<:Number,X} =
#            $poly{promote_type(T, S),X}

        $poly(coeffs::AbstractVector{T}, var::SymbolLike=Var(:x)) where {T} =
            $poly{T, Symbol(var)}(coeffs)
        $poly{T}(x::AbstractVector{S}, var::SymbolLike=Var(:x)) where {T,S} =
            $poly{T,Symbol(var)}(T.(x))

        function $poly{T}(coeffs::G, var::SymbolLike=Var(:x)) where {T,G}
            !Base.isiterable(G) && throw(ArgumentError("coeffs is not iterable"))
            cs = collect(T, coeffs)
            $poly{T, Symbol(var)}(cs)
        end
        function $poly(coeffs::G, var::SymbolLike=Var(:x)) where {G}
            !Base.isiterable(G) && throw(ArgumentError("coeffs is not iterable"))
            cs = collect(promote(coeffs...))
            $poly{eltype(cs), Symbol(var)}(cs)
        end

        $poly{T,X}(c::AbstractPolynomial{S,Y}) where {T,X,S,Y} = convert($poly{T,X}, c)
        $poly{T}(c::AbstractPolynomial{S,Y}) where {T,S,Y} = convert($poly{T}, c)
        $poly(c::AbstractPolynomial{S,Y}) where {S,Y} = convert($poly, c)

        $poly{T,X}(n::S) where {T, X, S<:Number} =
            T(n) *  one($poly{T, X})
        $poly{T}(n::S, var::SymbolLike = Var(:x)) where {T, S<:Number} =
            T(n) *  one($poly{T, Symbol(var)})
        $poly(n::S, var::SymbolLike = Var(:x))  where {S  <: Number} = n * one($poly{S, Symbol(var)})

        $poly{T}(var::SymbolLike=Var(:x)) where {T} = variable($poly{T, Symbol(var)})
        $poly(var::SymbolLike=Var(:x)) = variable($poly, Symbol(var))

        (p::$poly)(x::AbstractPolynomial) = polynomial_composition(p, x)
        (p::$poly)(x) = evalpoly(x, p)
    end
end
