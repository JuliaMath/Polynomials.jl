export ChebyshevT

"""
    ChebyshevT{<:Number}(coeffs::AbstractVector, var=:x)

Chebyshev polynomial of the first kind.

Construct a polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

# Examples

```jldoctest ChebyshevT
julia> using Polynomials

julia> ChebyshevT([1, 0, 3, 4])
ChebyshevT(1⋅T_0(x) + 3⋅T_2(x) + 4⋅T_3(x))

julia> ChebyshevT([1, 2, 3, 0], :s)
ChebyshevT(1⋅T_0(s) + 2⋅T_1(s) + 3⋅T_2(s))

julia> one(ChebyshevT)
ChebyshevT(1.0⋅T_0(x))
```
"""
struct ChebyshevT{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function ChebyshevT{T}(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

@register ChebyshevT

function ChebyshevT{T}(coeffs::OffsetArray{T,1, Array{T, 1}}, var::SymbolLike=:x) where {T <: Number}
    cs = zeros(T, 1 + lastindex(coeffs))
    cs[1 .+ (firstindex(coeffs):lastindex(coeffs))] = coeffs.parent
    ChebyshevT{T}(cs, var)
end


function Base.convert(P::Type{<:Polynomial}, ch::ChebyshevT)
    if length(ch) < 3
        return P(ch.coeffs, ch.var)
    end
    c0 = P(ch[end - 1], ch.var)
    c1 = P(ch[end], ch.var)
    @inbounds for i in degree(ch):-1:2
        tmp = c0
        c0 = P(ch[i - 2], ch.var) - c1
        c1 = tmp + c1 * variable(P) * 2
    end
    return c0 + c1 * variable(P)
end

function Base.convert(C::Type{<:ChebyshevT}, p::Polynomial)
    res = zero(C)
    @inbounds for i in degree(p):-1:0
        res = variable(C) * res + p[i]
    end
    return res
end


domain(::Type{<:ChebyshevT}) = Interval(-1, 1)
variable(P::Type{<:ChebyshevT}, var::SymbolLike=:x ) = P([0,1], var)
"""
    (::ChebyshevT)(x)

Evaluate the Chebyshev polynomial at `x`. If `x` is outside of the domain of [-1, 1], an error will be thrown. The evaluation uses Clenshaw Recursion.

# Examples
```jldoctest ChebyshevT
julia> using Polynomials

julia> c = ChebyshevT([2.5, 1.5, 1.0])
ChebyshevT(2.5⋅T_0(x) + 1.5⋅T_1(x) + 1.0⋅T_2(x))

julia> c(0)
1.5

julia> c.(-1:0.5:1)
5-element Array{Float64,1}:
 2.0
 1.25
 1.5
 2.75
 5.0
```
"""
function (ch::ChebyshevT{T})(x::S) where {T,S}
    x ∉ domain(ch) && error("$x outside of domain")
    R = promote_type(T, S)
    length(ch) == 0 && return zero(R)
    length(ch) == 1 && return R(ch[0])
    c0 = ch[end - 1]
    c1 = ch[end]
    @inbounds for i in lastindex(ch) - 2:-1:0
        c0, c1 = ch[i] - c1, c0 + c1 * 2x
    end
    return R(c0 + c1 * x)
end

function vander(P::Type{<:ChebyshevT}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    if n > 0
        A[:, 2] .= x
        @inbounds for i in 3:n + 1
            A[:, i] .= A[:, i - 1] .* 2x .- A[:, i - 2]
        end
    end
    return A
end

function integrate(p::ChebyshevT{T}, C::S) where {T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(C)
        return ChebyshevT([NaN])
    end
    n = length(p)
    if n == 1
        return ChebyshevT{R}([C, p[0]])
    end
    a2 = Vector{R}(undef, n + 1)
    a2[1] = zero(R)
    a2[2] = p[0]
    a2[3] = p[1] / 4
    @inbounds for i in 2:n - 1
        a2[i + 2] = p[i] / (2 * (i + 1))
        a2[i] -= p[i] / (2 * (i - 1))
    end
    a2[1] += R(C) - ChebyshevT(a2)(0)
    return ChebyshevT(a2, p.var)
end


function derivative(p::ChebyshevT{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    R  = eltype(one(T)/1)
    order == 0 && return convert(ChebyshevT{R}, p)
    hasnan(p) && return ChebyshevT(R[NaN], p.var)
    order > length(p) && return zero(ChebyshevT{R})


    q =  convert(ChebyshevT{R}, copy(p))
    n = length(p)
    der = Vector{R}(undef, n)

    for j in n:-1:2
        der[j] = 2j * q[j]
        q[j - 2] += j * q[j] / (j - 2)
    end
    if n > 1
        der[2] = 4q[2]
    end
    der[1] = q[1]

    pp = ChebyshevT(der, p.var)
    return order > 1 ?  derivative(pp, order - 1) : pp

end

function companion(p::ChebyshevT{T}) where T
    d = length(p) - 1
    d < 1 && error("Series must have degree greater than 1")
    d == 1 && return diagm(0 => [-p[0] / p[1]])
    R = eltype(one(T) / one(T))

    scl = vcat(1.0, fill(R(√0.5), d - 1))

    diag = vcat(√0.5, fill(R(0.5), d - 2))
    comp = diagm(1 => diag,
                 -1 => diag)
    monics = p.coeffs ./ p.coeffs[end]
    comp[:, end] .-= monics[1:d] .* scl ./ scl[end] ./ 2
    return R.(comp)
end

function Base.:+(p1::ChebyshevT, p2::ChebyshevT)
    p1.var != p2.var && error("Polynomials must have same variable")
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n]
    return ChebyshevT(c, p1.var)
end


function Base.:*(p1::ChebyshevT{T}, p2::ChebyshevT{S}) where {T,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    z1 = _c_to_z(p1.coeffs)
    z2 = _c_to_z(p2.coeffs)
    prod = fastconv(z1, z2)
    ret = ChebyshevT(_z_to_c(prod), p1.var)
    return truncate!(ret)
end

function Base.divrem(num::ChebyshevT{T}, den::ChebyshevT{S}) where {T,S}
    num.var != den.var && error("Polynomials must have same variable")
    n = length(num) - 1
    m = length(den) - 1

    R = typeof(one(T) / one(S))
    P = ChebyshevT{R}

    if n < m
        return zero(P), convert(P, num)
    elseif m == 0
        den[0] ≈ 0 && throw(DivideError())
        return num ./ den[end], zero(P)
    end

    znum = _c_to_z(num.coeffs)
    zden = _c_to_z(den.coeffs)
    quo, rem = _z_division(znum, zden)
    q_coeff = _z_to_c(quo)
    r_coeff = _z_to_c(rem)
    return P(q_coeff, num.var), P(r_coeff, num.var)
end

function showterm(io::IO, ::Type{ChebyshevT{T}}, pj::T, var, j, first::Bool, mimetype) where {N, T}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, hasneg(T)  && isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅T_$j($var)")
    return true
end


#=
zseries =#

function _c_to_z(cs::AbstractVector{T}) where {T}
    n = length(cs)
    U = typeof(one(T) / 2)
    zs = zeros(U, 2n - 1)
    zs[n:end] = cs ./ 2
    return zs .+ reverse(zs)
end

function _z_to_c(z::AbstractVector{T}) where {T}
    n = (length(z) + 1) ÷ 2
    cs = z[n:end]
    cs[2:n] *= 2
    return cs
end

function _z_division(z1::AbstractVector{T}, z2::AbstractVector{S}) where {T,S}
    R = eltype(one(T) / one(S))
    length(z1)
    length(z2)
    if length(z2) == 1
        z1 ./= z2
        return z1, zero(R)
    elseif length(z1) < length(z2)
        return zero(R), R.(z1)
    end
    dlen = length(z1) - length(z2)
    scl = z2[1]
    z2 ./= scl
    quo = Vector{R}(undef, dlen + 1)
    i = 1
    j = dlen + 1
    while i < j
        r = z1[i]
        quo[i] = z1[i]
        quo[end - i + 1] = r
        tmp = r .* z2
        z1[i:i + length(z2) - 1] .-= tmp
        z1[j:j + length(z2) - 1] .-= tmp
        i += 1
        j -= 1
    end

    r = z1[i]
    quo[i] = r
    tmp = r * z2
    z1[i:i + length(z2) - 1] .-= tmp
    quo ./= scl
    rem = z1[i + 1:i - 2 + length(z2)]
    return quo, rem
end
