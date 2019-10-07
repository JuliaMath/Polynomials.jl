export AbstractPolynomial

const SymbolLike = Union{AbstractString,Char,Symbol}

"""
    AbstractPolynomial{<:Number}

An abstract container for various polynomials. 

# Properties
- `coeffs` - The coefficients of the polynomial
- `var` - The indeterminate of the polynomial
"""
abstract type AbstractPolynomial{T <: Number} end


macro register(name)
    poly = esc(name)
    quote
    Base.convert(::Type{P}, p::P) where {P <: $poly} = p
    Base.convert(P::Type{<:$poly}, p::$poly) where {T} = P(p.coeffs, p.var)
    Base.promote_rule(::Type{$poly{T}}, ::Type{$poly{S}}) where {T,S} = $poly{promote_type(T, S)}
    Base.promote_rule(::Type{$poly{T}}, ::Type{S}) where {T,S <: Number} = $poly{promote_type(T, S)}

    function (p::$poly)(x::AbstractVector)
        Base.depwarn("Calling p(x::AbstractVector is deprecated. Use p.(x) instead.", Symbol("(p::AbstractPolynomial)"))
        return p.(x)
    end

    $poly(coeffs::AbstractVector{T}, var::SymbolLike = :x) where {T} = $poly{T}(coeffs, Symbol(var))
    $poly(n::Number, var = :x) = $poly([n], var)
    $poly{T}(n::S, var = :x) where {T,S <: Number} = $poly(T(n), var)
    $poly{T}(x::AbstractVector{S}, var = :x) where {T,S <: Number} = $poly(T.(x), var)
    $poly(x, var = :x) = $poly(collect(x), var)

    end
end
