export AbstractPolynomial



const SymbolLike = Union{AbstractString,Char,Symbol}

"""
    AbstractPolynomial{T,X}

An abstract type for various polynomials.

A polynomial type holds an indeterminate `X`; coefficients of type `T`, stored in some container type; and an implicit basis, such as the standard basis.

# Properties
- `coeffs` - The coefficients of the polynomial

# The type `T`

`T` need not be `T <: Number`, at the moment it is not restricted

Some `T`s will not be successful

* scalar mult: `c::Number * p::Polynomial` should be defined
* scalar mult: `c::T * p:Polynomial{T}` An  ambiguity when `T <: AbstractPolynomial`
* scalar mult: `p:Polynomial{T} * c::T` need not commute

* scalar add/sub: `p::Polynomial{T} + q::Polynomial{T}` should be defined
* scalar sub: `p::Polynomial{T} - p::Polynomial{T}`  generally needs `zeros(T,1)` defined for `zero(Polynomial{T})`

* poly mult: `p::Polynomial{T} * q::Polynomial{T}` Needs "`T * T`" defined (e.g. `Base.promote_op(*, Vector{Int}, Vector{Int}))` needs to be something.)
* poly powers: `p::Polynomial{T}^2` needs "`T^2`" defined

* implicit promotion: `p::Polynomial{T} + c::Number`  needs `convert(T, c)` defined
* implicit promotion: `p::Polynomial{T} + c::T`  ambiguity if `T <: AbstractPolynomial`

* evaluation: `p::Polynomial{T}(s::Number)`
* evaluation `p::Polynomial{T}(c::T)`   needs `T*T` defined

* derivatives: `derivative(p::Polynomial{T})`
* integrals: `integrate(p::Polynomial{T})`


"""
abstract type AbstractPolynomial{T,X} end


## -----
# We want  ⟒(P{α…,T}) = P{α…}; this default
# works for most cases
⟒(P::Type{<:AbstractPolynomial}) = constructorof(P)

# convert `as` into polynomial of type P based on instance, inheriting variable
# (and for LaurentPolynomial the offset)
_convert(p::P, as) where {T,X,P <: AbstractPolynomial{T,X}} = ⟒(P)(as, X)


"""
    Polynomials.@register(name)

Given a polynomial with `name`, creates some common convenience constructors and conversions to minimize code required for implementation of a new polynomial type.

# Example
```julia
struct MyPolynomial{T,X} <: AbstractPolynomial{T,X} end

Polynomials.@register MyPolynomial
```

# Implementations
This will implement simple self-conversions like `convert(::Type{MyPoly}, p::MyPoly) = p` and creates two promote rules. The first allows promotion between two types (e.g. `promote(Polynomial, ChebyshevT)`) and the second allows promotion between parametrized types (e.g. `promote(Polynomial{T}, Polynomial{S})`).

For constructors, it implements the shortcut for `MyPoly(...) = MyPoly{T}(...)`, singleton constructor `MyPoly(x::Number, ...)`,  conversion constructor `MyPoly{T}(n::S, ...)`, and `variable` alternative  `MyPoly(var=:x)`.
"""
macro register(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        function Base.convert(P::Type{<:$poly}, p::$poly{T,X}) where {T,X}
            isconstant(p) && return constructorof(P){eltype(P),indeterminate(P)}(constantterm(p))
            constructorof(P){eltype(P), indeterminate(P,p)}(_coeffs(p))
        end
        Base.promote(p::P, q::Q) where {X, T, P <:$poly{T,X}, Q <: $poly{T,X}} = p,q
        Base.promote_rule(::Type{<:$poly{T,X}}, ::Type{<:$poly{S,X}}) where {T,S,X} =  $poly{promote_type(T, S),X}
        Base.promote_rule(::Type{<:$poly{T,X}}, ::Type{S}) where {T,S<:Number,X} =
            $poly{promote_type(T, S),X}
        $poly(coeffs::AbstractVector{T}, var::SymbolLike = :x) where {T} =
            $poly{T, Symbol(var)}(coeffs)
        $poly{T}(x::AbstractVector{S}, var::SymbolLike = :x) where {T,S} =
            $poly{T,Symbol(var)}(T.(x))
        function $poly(coeffs::G, var::SymbolLike=:x) where {G}
            !Base.isiterable(G) && throw(ArgumentError("coeffs is not iterable"))
            cs = collect(coeffs)
            $poly{eltype(cs), Symbol(var)}(cs)
        end
        $poly{T,X}(c::AbstractPolynomial{S,Y}) where {T,X,S,Y} = convert($poly{T,X}, c)
        $poly{T}(c::AbstractPolynomial{S,Y}) where {T,S,Y} = convert($poly{T}, c)
        $poly(c::AbstractPolynomial{S,Y}) where {S,Y} = convert($poly, c)
        $poly{T,X}(n::S) where {T, X, S<:Number} =
            T(n) *  one($poly{T, X})
        $poly{T}(n::S, var::SymbolLike = :x) where {T, S<:Number} =
            T(n) *  one($poly{T, Symbol(var)})
        $poly(n::S, var::SymbolLike = :x)  where {S  <: Number} = n * one($poly{S, Symbol(var)})
        $poly{T}(var::SymbolLike=:x) where {T} = variable($poly{T, Symbol(var)})
        $poly(var::SymbolLike=:x) = variable($poly, Symbol(var))
        (p::$poly)(x) = evalpoly(x, p)
    end
end


macro registerN(name, params...)
    poly = esc(name)
    αs = tuple(esc.(params)...)
    quote
        Base.convert(::Type{P}, q::Q) where {$(αs...),T, P<:$poly{$(αs...),T}, Q <: $poly{$(αs...),T}} = q
        Base.convert(::Type{$poly{$(αs...)}}, q::Q) where {$(αs...),T, Q <: $poly{$(αs...),T}} = q
        Base.promote(p::P, q::Q) where {$(αs...),T, X, P <:$poly{$(αs...),T,X}, Q <: $poly{$(αs...),T,X}} = p,q
        Base.promote_rule(::Type{<:$poly{$(αs...),T,X}}, ::Type{<:$poly{$(αs...),S,X}}) where {$(αs...),T,S,X} =
            $poly{$(αs...),promote_type(T, S),X}
        Base.promote_rule(::Type{<:$poly{$(αs...),T,X}}, ::Type{S}) where {$(αs...),T,X,S<:Number} =
            $poly{$(αs...),promote_type(T,S),X}

        function $poly{$(αs...),T}(x::AbstractVector{S}, var::SymbolLike = :x) where {$(αs...),T,S}
            $poly{$(αs...),T, Symbol(var)}(T.(x))
        end
        $poly{$(αs...)}(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {$(αs...),T} =
            $poly{$(αs...),T,Symbol(var)}(coeffs)
        $poly{$(αs...),T,X}(c::AbstractPolynomial{S,Y}) where {$(αs...),T,X,S,Y} = convert($poly{$(αs...),T,X}, c)
        $poly{$(αs...),T}(c::AbstractPolynomial{S,Y}) where {$(αs...),T,S,Y} = convert($poly{$(αs...),T}, c)
        $poly{$(αs...),}(c::AbstractPolynomial{S,Y}) where {$(αs...),S,Y} = convert($poly{$(αs...),}, c)

        $poly{$(αs...),T,X}(n::Number) where {$(αs...),T,X} = T(n)*one($poly{$(αs...),T,X})
        $poly{$(αs...),T}(n::Number, var::SymbolLike = :x) where {$(αs...),T} = T(n)*one($poly{$(αs...),T,Symbol(var)})
        $poly{$(αs...)}(n::S, var::SymbolLike = :x) where {$(αs...), S<:Number} =
            n*one($poly{$(αs...),S,Symbol(var)})
        $poly{$(αs...),T}(var::SymbolLike=:x) where {$(αs...), T} = variable($poly{$(αs...),T,Symbol(var)})
        $poly{$(αs...)}(var::SymbolLike=:x) where {$(αs...)} = variable($poly{$(αs...)},Symbol(var))
        (p::$poly)(x) = evalpoly(x, p)
    end
end
