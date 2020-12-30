export Polynomial

"""
    Polynomial{T<:Number}(coeffs::AbstractVector{T}, [var = :x])

Construct a polynomial from its coefficients `coeffs`, lowest order first, optionally in
terms of the given variable `var` which may be a character, symbol, or a string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct this through
`Polynomial([a_0, a_1, ..., a_n])`.

The usual arithmetic operators are overloaded to work with polynomials as well as
with combinations of polynomials and scalars. However, operations involving two
polynomials of different variables causes an error except those involving a constant polynomial.

!!! note
    `Polynomial` is not axis-aware, and it treats `coeffs` simply as a list of coefficients with the first 
    index always corresponding to the constant term. In order to use the axis of `coeffs` as exponents, 
    consider using a [`LaurentPolynomial`](@ref) or possibly a [`SparsePolynomial`](@ref).

# Examples
```jldoctest
julia> Polynomial([1, 0, 3, 4])
Polynomial(1 + 3*x^2 + 4*x^3)

julia> Polynomial([1, 2, 3], :s)
Polynomial(1 + 2*s + 3*s^2)

julia> one(Polynomial)
Polynomial(1.0)
```
"""
struct Polynomial{T <: Number} <: StandardBasisPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Polynomial{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        if Base.has_offset_axes(coeffs)
          @warn "ignoring the axis offset of the coefficient vector"
        end
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        c = OffsetArrays.no_offset_view(coeffs) # ensure 1-based indexing
        last_nz = findlast(!iszero, c)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(c[1:last], var)
    end
end

@register Polynomial

"""
    (p::Polynomial)(x)

Evaluate the polynomial using [Horner's Method](https://en.wikipedia.org/wiki/Horner%27s_method), also known as synthetic division, as implemented in `evalpoly` of base `Julia`.

# Examples
```jldoctest
julia> p = Polynomial([1, 0, 3])
Polynomial(1 + 3*x^2)

julia> p(0)
1

julia> p.(0:3)
4-element Array{Int64,1}:
  1
  4
 13
 28
```
"""
(p::Polynomial{T})(x::S) where {T,S} = evalpoly(x, coeffs(p))


   
function Base.:+(p1::Polynomial{T}, p2::Polynomial{S}) where {T, S}
    n1, n2 = length(p1), length(p2)
    if n1 > 1 && n2 > 1
       p1.var != p2.var && error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    c = zeros(R, max(n1, n2))
    if n1 > 1 && n2 > 1
       if n1 >= n2
          c .= p1.coeffs
          for i = eachindex(p2.coeffs)
            c[i] += p2.coeffs[i]
          end
        else
            c .= p2.coeffs
            for i = eachindex(p1.coeffs)
              c[i] += p1.coeffs[i]
            end
        end
        return Polynomial(c, p1.var)
    elseif n1 <= 1
      c .= p2.coeffs
      c[1] += p1[0]
      return Polynomial(c, p2.var)
    else 
      c .= p1.coeffs
      c[1] += p2[0]
      return Polynomial(c, p1.var)
    end
end


function Base.:*(p1::Polynomial{T}, p2::Polynomial{S}) where {T,S}

    n, m = length(p1)-1, length(p2)-1 # not degree, so pNULL works
    if n > 0 && m > 0
        p1.var != p2.var && error("Polynomials must have same variable")
        R = promote_type(T, S)
        c = zeros(R, m + n + 1)
        for i in 0:n, j in 0:m
            @inbounds c[i + j + 1] += p1[i] * p2[j]
        end
        return Polynomial(c, p1.var)
    elseif n <= 0
        return Polynomial(p2.coeffs * p1[0], p2.var)
    else
        return Polynomial(p1.coeffs * p2[0], p1.var)
    end

end
