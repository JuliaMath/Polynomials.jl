export AbstractPolynomial

const SymbolLike = Union{AbstractString,Char,Symbol}

"""
    AbstractPolynomial{<:Number}

An abstract container for various polynomials. 

# Properties
- `coeffs` - The coefficients of the polynomial
- `var` - The indeterminate of the polynomial
"""
abstract type AbstractPolynomial{T<:Number} end


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

For constructors, it implements the shortcut for `MyPoly(...) = MyPoly{T}(...)`, singleton constructor `MyPoly(x::Number, ...)`, and conversion constructor `MyPoly{T}(n::S, ...)`.
"""
macro register(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.convert(P::Type{<:$poly}, p::$poly) where {T} = P(coeffs(p), p.var)
        Base.promote_rule(::Type{$poly{T}}, ::Type{$poly{S}}) where {T,S} =
            $poly{promote_type(T, S)}
        Base.promote_rule(::Type{$poly{T}}, ::Type{S}) where {T,S<:Number} =
            $poly{promote_type(T, S)}

        function (p::$poly)(x::AbstractVector)
            throw(ArgumentError("Calling p(x::AbstractVector is not supported. Use p.(x) instead."))
        end

        $poly(coeffs::AbstractVector{T}, var::SymbolLike = :x) where {T} =
            $poly{T}(coeffs, Symbol(var))
        $poly(n::Number, var = :x) = $poly([n], var)
        $poly{T}(n::S, var = :x) where {T,S<:Number} = $poly(T(n), var)
        $poly{T}(x::AbstractVector{S}, var = :x) where {T,S<:Number} = $poly(T.(x), var)
    end
end
