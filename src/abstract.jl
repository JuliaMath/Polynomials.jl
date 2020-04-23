export AbstractPolynomial

const SymbolLike = Union{AbstractString,Char,Symbol}

"""
    AbstractPolynomial{T}

An abstract container for various polynomials. 

# Properties
- `coeffs` - The coefficients of the polynomial
- `var` - The indeterminate of the polynomial
"""
abstract type AbstractPolynomial{T} end


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
        Base.convert(P::Type{<:$poly}, p::$poly{T}) where {T} = P(coeffs(p), p.var)
        Base.promote_rule(::Type{$poly{T}}, ::Type{$poly{S}}) where {T,S} =
            $poly{promote_type(T, S)}
        Base.promote_rule(::Type{$poly{T}}, ::Type{S}) where {T,S<:Number} =
            $poly{promote_type(T, S)}
        $poly(coeffs::AbstractVector{T}, var::SymbolLike = :x) where {T} =
            $poly{T}(coeffs, Symbol(var))
        $poly{T}(x::AbstractVector{S}, var = :x) where {T,S<:Number} =
            $poly(T.(x), var)
        $poly{T}(n::S, var = :x) where {T, S<:Number} =
            $poly(T[n], var)
        $poly(n::Number, var = :x) = $poly([n], var)
        $poly{T}(var::SymbolLike=:x) where {T} = variable($poly{T}, Symbol(var))
        $poly(var::SymbolLike=:x) = variable($poly, Symbol(var))
    end
end


# Macros to register POLY{α, T} and POLY{α, β, T}
macro register1(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote_rule(::Type{$poly{α,T}}, ::Type{$poly{α,S}}) where {α,T,S} =
            $poly{α,promote_type(T, S)}
        Base.promote_rule(::Type{$poly{α,T}}, ::Type{S}) where {α,T,S<:Number} = 
            $poly{α,promote_type(T,S)}
        function $poly{α,T}(x::AbstractVector{S}, var::Polynomials.SymbolLike = :x) where {α,T,S}
            $poly{α,T}(T.(x), Symbol(var))
        end
        $poly{α}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike=:x) where {α,T} =
            $poly{α,T}(coeffs, Symbol(var))
        $poly{α,T}(n::Number, var::Polynomials.SymbolLike = :x) where {α,T} = n*one($poly{α,T}, Symbol(var))
        $poly{α}(n::Number, var::Polynomials.SymbolLike = :x) where {α} = n*one($poly{α}, Symbol(var))
        $poly{α,T}(var::Polynomials.SymbolLike=:x) where {α, T} = variable($poly{α,T}, Symbol(var))
        $poly{α}(var::Polynomials.SymbolLike=:x) where {α} = variable($poly{α}, Symbol(var))
    end
end


# Macro to register POLY{α, β, T}
macro register2(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote_rule(::Type{$poly{α,β,T}}, ::Type{$poly{α,β,S}}) where {α,β,T,S} =
            $poly{α,β,promote_type(T, S)}
        Base.promote_rule(::Type{$poly{α,β,T}}, ::Type{S}) where {α,β,T,S<:Number} =
            $poly{α,β,promote_type(T, S)}
        $poly{α,β}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where {α,β,T} =
            $poly{α,β,T}(coeffs, Symbol(var))
        $poly{α,β,T}(x::AbstractVector{S}, var = :x) where {α,β,T,S<:Number} = $poly{α,β,T}(T.(x), var)
        $poly{α,β,T}(n::Number, var = :x) where {α,β,T} = n*one($poly{α,β,T}, var)
        $poly{α,β}(n::Number, var = :x) where {α,β} = n*one($poly{α,β}, var)
        $poly{α,β,T}(var=:x) where {α,β, T} = variable($poly{α,β,T}, var)
        $poly{α,β}(var=:x) where {α,β} = variable($poly{α,β}, var)
    end
end

