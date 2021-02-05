export AbstractPolynomial

const SymbolLike = Union{AbstractString,Char,Symbol}

"""
    AbstractPolynomial{T, X}

An abstract container for various polynomials. 

# Properties
- `coeffs` - The coefficients of the polynomial
"""
abstract type AbstractPolynomial{T, X} end

# We want  ⟒(P{α…,T}) = P{α…}; this default
# works for most cases
⟒(P::Type{<:AbstractPolynomial}) = constructorof(P)

# convert `as` into polynomial of type P based on instance, inheriting variable
# (and for LaurentPolynomial the offset)
_convert(p::P, as) where {P <: AbstractPolynomial} = ⟒(P)(as, var(P))

"""
    Polynomials.@register(name)

Given a polynomial with `name`, creates some common convenience constructors and conversions to minimize code required for implementation of a new polynomial type.

# Example
```julia
struct MyPolynomial{T} <: AbstractPolynomial{T} end

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
        Base.convert(P::Type{<:$poly}, p::$poly{T}) where {T} = P(coeffs(p), var(p))
        Base.promote(p::P, q::Q) where {X, T, P <:$poly{T,X}, Q <: $poly{T,X}} = p,q
        Base.promote_rule(::Type{<:$poly{T,X}}, ::Type{<:$poly{S,X}}) where {T,S,X} =
            $poly{promote_type(T, S,X)}
        Base.promote_rule(::Type{<:$poly{T,X}}, ::Type{S}) where {T,S<:Number,X} =
            $poly{promote_type(T, S),X}
        $poly(coeffs::AbstractVector{T}, var::SymbolLike = :x) where {T} =
            $poly{T, Symbol(var)}(coeffs)
        $poly{T}(x::AbstractVector{S}, var::SymbolLike = :x) where {T,S<:Number} =
            $poly(T.(x), Symbol(var))
        function $poly(coeffs::G, var::SymbolLike=:x) where {G}
            !Base.isiterable(G) && throw(ArgumentError("coeffs is not iterable"))
            $poly(collect(coeffs), var)
        end
        $poly{T,X}(n::S) where {X, T, S<:Number} =
            n *  one($poly{T}, X)
        $poly{T}(n::S, var::SymbolLike = :x) where {T, S<:Number} =
            n *  one($poly{T}, Symbol(var))
        $poly(n::S, var::SymbolLike = :x)  where {S  <: Number} = n * one($poly{S}, Symbol(var))
        $poly{T}(var::SymbolLike=:x) where {T} = variable($poly{T}, Symbol(var))
        $poly(var::SymbolLike=:x) = variable($poly, Symbol(var))
    end
end


macro registerN(name, params...)
    poly = esc(name)
    αs = tuple(esc.(params)...)
    quote
        Base.convert(::Type{P}, q::Q) where {$(αs...),T, P<:$poly{$(αs...),T}, Q <: $poly{$(αs...),T}} = q
        Base.convert(::Type{$poly{$(αs...)}}, q::Q) where {$(αs...),T, Q <: $poly{$(αs...),T}} = q        
        Base.promote(p::P, q::Q) where {$(αs...),T, P <:$poly{$(αs...),T}, Q <: $poly{$(αs...),T}} = p,q
        Base.promote_rule(::Type{<:$poly{$(αs...),T}}, ::Type{<:$poly{$(αs...),S}}) where {$(αs...),T,S} =
            $poly{$(αs...),promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{$(αs...),T}}, ::Type{S}) where {$(αs...),T,S<:Number} = 
            $poly{$(αs...),promote_type(T,S)}

        function $poly{$(αs...),T}(x::AbstractVector{S}, var::SymbolLike = :x) where {$(αs...),T,S}
            $poly{$(αs...),T}(T.(x), Symbol(var))
        end
        $poly{$(αs...)}(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {$(αs...),T} =
            $poly{$(αs...),T}(coeffs, Symbol(var))
        $poly{$(αs...),T}(n::Number, var::SymbolLike = :x) where {$(αs...),T} = n*one($poly{$(αs...),T}, Symbol(var))
        $poly{$(αs...)}(n::Number, var::SymbolLike = :x) where {$(αs...)} = n*one($poly{$(αs...)}, Symbol(var))
        $poly{$(αs...),T}(var::SymbolLike=:x) where {$(αs...), T} = variable($poly{$(αs...),T}, Symbol(var))
        $poly{$(αs...)}(var::SymbolLike=:x) where {$(αs...)} = variable($poly{$(αs...)}, Symbol(var))
    end
end


# deprecated. If desired,  replace with  @registerN  type  parameters... macro
# Macros to register POLY{α, T} and POLY{α, β, T}
macro register1(name)
    @warn "@register1 is deprecated use @registerN"
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote(p::P, q::Q) where {α,T, P <:$poly{α,T}, Q <: $poly{α,T}} = p,q
        Base.promote_rule(::Type{<:$poly{α,T}}, ::Type{<:$poly{α,S}}) where {α,T,S} =
            $poly{α,promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{α,T}}, ::Type{S}) where {α,T,S<:Number} = 
            $poly{α,promote_type(T,S)}
        function $poly{α,T}(x::AbstractVector{S}, var::SymbolLike = :x) where {α,T,S}
            $poly{α,T}(T.(x), Symbol(var))
        end
        $poly{α}(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {α,T} =
            $poly{α,T}(coeffs, Symbol(var))
        $poly{α,T}(n::Number, var::SymbolLike = :x) where {α,T} = n*one($poly{α,T}, Symbol(var))
        $poly{α}(n::Number, var::SymbolLike = :x) where {α} = n*one($poly{α}, Symbol(var))
        $poly{α,T}(var::SymbolLike=:x) where {α, T} = variable($poly{α,T}, Symbol(var))
        $poly{α}(var::SymbolLike=:x) where {α} = variable($poly{α}, Symbol(var))
    end
end


# Macro to register POLY{α, β, T}
macro register2(name)
    @warn "@register2  is deprecated use @registerN"
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote(p::P, q::Q) where {α,β,T, P <:$poly{α,β,T}, Q <: $poly{α,β,T}} = p,q        
        Base.promote_rule(::Type{<:$poly{α,β,T}}, ::Type{<:$poly{α,β,S}}) where {α,β,T,S} =
            $poly{α,β,promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{α,β,T}}, ::Type{S}) where {α,β,T,S<:Number} =
            $poly{α,β,promote_type(T, S)}
        $poly{α,β}(coeffs::AbstractVector{T}, var::SymbolLike = :x) where {α,β,T} =
            $poly{α,β,T}(coeffs, Symbol(var))
        $poly{α,β,T}(x::AbstractVector{S}, var::SymbolLike = :x) where {α,β,T,S<:Number} = $poly{α,β,T}(T.(x), var)
        $poly{α,β,T}(n::Number, var::SymbolLike = :x) where {α,β,T} = n*one($poly{α,β,T}, Symbol(var))
        $poly{α,β}(n::Number, va::SymbolLiker = :x) where {α,β} = n*one($poly{α,β}, Symbol(var))
        $poly{α,β,T}(var::SymbolLike=:x) where {α,β, T} = variable($poly{α,β,T}, Symbol(var))
        $poly{α,β}(var::SymbolLike=:x) where {α,β} = variable($poly{α,β}, Symbol(var))
    end
end

