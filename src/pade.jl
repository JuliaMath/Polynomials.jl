module PadeApproximation

using ..Polynomials
export Pade, padeval

#=
Pade approximation

Note: This can be moved to polynomials/Polynomial.jl after Poly type is removed
=#

export Pade


"""
    Pade(::Polynomial, m::Integer, n::Integer)
    Pade(::Polynomial, ::Polynomial)

Return Pade approximation of polynomial.

# References
[Pade Approximant](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant)
"""
struct Pade{T <: Number,S <: Number}
    p::Union{Poly{T}, Polynomial{T}}
    q::Union{Poly{S}, Polynomial{S}}
    var::Symbol
    function Pade{T,S}(p::Union{Poly{T}, Polynomial{T}}, q::Union{Poly{S}, Polynomial{S}}) where {T,S}
        if p.var != q.var error("Polynomials must have same variable") end
        new{T,S}(p, q, p.var)
    end
end

Pade(p::Polynomial{T}, q::Polynomial{S}) where {T <: Number,S <: Number} = Pade{T,S}(p, q)

function Pade(c::Polynomial{T}, m::Integer, n::Integer) where {T}
    m + n < length(c) || error("m + n must be less than the length of the Polynomial")
    rold = Polynomial([zeros(T, m + n + 1);one(T)], c.var)
    rnew = Polynomial(c[0:m + n], c.var)
    uold = Polynomial([one(T)], c.var)
    vold = Polynomial([zero(T)], c.var)
    unew, vnew = vold, uold
    @inbounds for i = 1:n
        temp0, temp1, temp2 = rnew, unew, vnew
        q, rnew = divrem(rold, rnew)
        unew, vnew = uold - q * unew, vold - q * vnew
        rold, uold, vold = temp0, temp1, temp2

    end
    if vnew[0] == 0
        d = gcd(rnew, vnew)
        rnew = rnew ÷ d
        vnew = vnew ÷ d
    end
    Pade(rnew / vnew[0], vnew / vnew[0])
end


Pade(p::Poly{T}, q::Poly{S}) where {T <: Number,S <: Number} = Pade{T,S}(p, q)

function Pade(c::Poly{T}, m::Integer, n::Integer) where {T}
    m + n < length(c) || error("m + n must be less than the length of the polynomial")
    rold = Poly([zeros(T, m + n + 1);one(T)], c.var)
    rnew = Poly(c[0:m + n], c.var)
    uold = Poly([one(T)], c.var)
    vold = Poly([zero(T)], c.var)
    unew, vnew = vold, uold
    @inbounds for i = 1:n
        temp0, temp1, temp2 = rnew, unew, vnew
        q, rnew = divrem(rold, rnew)
        unew, vnew = uold - q * unew, vold - q * vnew
        rold, uold, vold = temp0, temp1, temp2

    end
    if vnew[0] == 0
        d = gcd(rnew, vnew)
        rnew = rnew ÷ d
        vnew = vnew ÷ d
    end
    Pade(rnew / vnew[0], vnew / vnew[0])
end

"""
    (::Pade)(x)

Evaluate the Pade approximant at the given point.

# Examples
```jldoctest
julia> using SpecialFunctions

julia> p = Polynomial(@.(1 // BigInt(gamma(1:17))));

julia> pade = Pade(p, 8, 8);

julia> pade(1.0) ≈ exp(1.0)
true

```
"""
(PQ::Pade)(x) = PQ.p(x) / PQ.q(x)


## Compat
padeval(PQ::Pade, x::Number) = PQ(x)
padeval(PQ::Pade, x) = PQ.(x)



end
