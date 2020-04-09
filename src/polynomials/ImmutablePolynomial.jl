export ImmutablePolynomial

"""
    ImmutablePolynomial{T<:Number}(coeffs::AbstractVector{T}, var=:x)

Construct an immutable polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct this through
`ImmutablePolynomial((a_0, a_1, ..., a_n))`.

The usual arithmetic operators are overloaded to work with polynomials as well as
with combinations of polynomials and scalars. However, operations involving two
polynomials of different variables causes an error.

This has the advantage over `Polynomial` as it can take advantage of faster polynomial evaluation
provided by `evalpoly` (borrowed from Julia 1.4).

# Examples

```jldoctest
julia> ImmutablePolynomial([1, 0, 3, 4])
ImmutablePolynomial(1 + 3*x^2 + 4*x^3)

julia> ImmutablePolynomial([1, 2, 3], :s)
ImmutablePolynomial(1 + 2*s + 3*s^2)

julia> one(ImmutablePolynomial)
ImmutablePolynomial(1.0)
```
"""
struct ImmutablePolynomial{N, T <: Number} <: AbstractPolynomial{T}
    coeffs::NTuple{N, T}
    var::Symbol
    function ImmutablePolynomial{N,T}(coeffs::NTuple{N,T}, var::Symbol=:x) where {N, T <: Number}
        new{N,T}(coeffs, var)
    end
    function ImmutablePolynomial{T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {T <: Number}
        length(coeffs) == 0 && return new{1, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        N = max(1, last_nz === nothing ? 0 : last_nz)
        return new{N,T}(NTuple{N, T}(x for x in coeffs[1:N]), var)
    end
   ImmutablePolynomial(coeffs::AbstractVector{T}, var::Symbol=:x) where {T <: Number} =  ImmutablePolynomial{T}(coeffs, var)
   function ImmutablePolynomial(coeffs::NTuple{N, T}, var::Symbol=:x) where {N, T}
        last_nz = findlast(!iszero, coeffs)
        M = last_nz === nothing ? 0 : last_nz
        if M < N
            new{M,T}(NTuple{M,T}(x  for x in coeffs[1:M]), var)
        else
            new{N,T}(coeffs, var)
        end
   end
    ImmutablePolynomial(x::Number, var::Symbol=:x) = ImmutablePolynomial((x,), var)
end

## Can't do this, as we keep the parameter N
## @register Polynomial
Base.convert(::Type{<:Polynomial}, p::ImmutablePolynomial{N,T}) where {N, T} = Polynomial{T}([coeffs(p)...], p.var)
Base.convert(::Type{ImmutablePolynomial{N,R}}, p::ImmutablePolynomial{N,S}) where {N, R, S} = ImmutablePolynomial(R.(coeffs(p)), p.var)
function Base.promote_rule(::Type{ImmutablePolynomial{N,T}},::Type{ImmutablePolynomial{M,S}}) where {N,T,M,S}
    NN = max(N,M)
    R = promote_type(T,S)
    ImmutablePolynomial{NN, R}
end
Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{S}) where {N,T,S} = ImmutablePolynomial{N, promote_type{T,S}}
ImmutablePolynomial{N, T}(x::S, var=:x) where  {N, T, S <: Number} = ImmutablePolynomial{T}(x, var)
ImmutablePolynomial{T}(x::S, var=:x) where  {T, S <: Number} = ImmutablePolynomial(T.(x), var)


                       
function showterm(io::IO, ::Type{ImmutablePolynomial{N,T}}, pj::T, var, j, first::Bool, mimetype) where {N,T}
    if iszero(pj) return false end
    pj = printsign(io, pj, first, mimetype)
    if !(pj == one(T) && !(showone(T) || j == 0))
        printcoefficient(io, pj, j, mimetype)
    end
    printproductsign(io, pj, j, mimetype)
    printexponent(io, var, j, mimetype)
    return true
end

domain(::Type{<:ImmutablePolynomial}) = Interval(-Inf, Inf)
mapdomain(::Type{<:ImmutablePolynomial}, x::AbstractArray) = x
degree(p::ImmutablePolynomial{N,T}) where {N, T} = N==1 && iszero(coeffs(p)[1]) ? -1 : N - 1

(p::ImmutablePolynomial{N, T})(x::S) where {N, T,S} = evalpoly(x, coeffs(p))

function fromroots(P::Type{<:ImmutablePolynomial}, r::AbstractVector{T}; var::SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j in 1:n, i in j:-1:1
        c[i + 1] = c[i + 1] - r[j] * c[i]
    end
    return ImmutablePolynomial(NTuple{n+1, T}(x for x in reverse(c)), var)
end

function vander(P::Type{<:ImmutablePolynomial}, x::AbstractVector{T}, n::Integer) where {T <: Number}
     A = Matrix{T}(undef, length(x), n + 1)
     A[:, 1] .= one(T)
     @inbounds for i in 1:n
         A[:, i + 1] = A[:, i] .* x
     end
     return A
end

function integrate(p::ImmutablePolynomial{N, T}, k::S) where {N, T,S <: Number}
    integrate(convert(Polynomial, p),  k)
end

function derivative(p::ImmutablePolynomial{T}, order::Integer = 1) where {T}
    order < 0 && error("Order of derivative must be non-negative")
    R = eltype(one(T)/1)
    order == 0 && return convert(Polynomial{R}, p)
    hasnan(p) && return Polynomial(R[NaN], p.var)
    order > length(p) && return zero(Polynomial{R})

    n = length(p)
    a2 = Vector{R}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return ImmutablePolynomial(NTuple{n-2, R}(a   for  a in a2), p.var)
end

companion(p::ImmutablePolynomial{T}) where T = companion(convert(Polynomial{T},p))
roots(p::ImmutablePolynomial{T}; kwargs...)  where  {T} = roots(convert(Polynomial{T}, p); kwargs...)

function Base.:+(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    R = promote_type(T,S)
    m = max(N, M)
    return ImmutablePolynomial{m, R}(NTuple{m, R}(p1[i] + p2[1] for i in 0:m-1), p1.var)
end


function Base.:+(p::ImmutablePolynomial{N, T}, c::S) where {N, T,S<:Number}
    R = promote_type(T, S)
    cs = R[pi for pi in coeffs(p)]
    cs[1] += c
    ImmutablePolynomial{N,R}(NTuple{N,R}(c for c in cs), p.var)
end

function Base.:*(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    p1.var != p2.var && error("Polynomials must have same variable")

    R = promote_type(T, S)
    cs = zeros(R, N + M -1)
    for i in 0:N-1, j in 0:M-1
        @inbounds cs[i + j + 1] += p1[i] * p2[j]
    end
    return ImmutablePolynomial{N+M-1,R}(NTuple{N+M-1,R}(c for c in cs), p1.var)
end

Base.:-(p::ImmutablePolynomial{N,T}) where {N,T} = ImmutablePolynomial(NTuple{N,T}(-pi for pi in p.coeffs), p.var)
Base.:*(c::Number, p::ImmutablePolynomial) = *(p, c)


function Base.:*(p::ImmutablePolynomial{N,T}, c::Number) where {N, T}
    R = eltype(one(T)*c)
    return ImmutablePolynomial(NTuple{N,R}(c*pi for pi in p.coeffs), p.var)
end


function Base.divrem(num::ImmutablePolynomial{N, T}, den::ImmutablePolynomial{M,S}) where {N,T,M,S}
    num.var != den.var && error("Polynomials must have same variable")
    n = degree(num)
    m = degree(den)
    m == -1 && throw(DivideError())
    R = typeof(one(T) / one(S))
    if m > n
        return ImmutablePolynomial{1,R}((0,)), convert(ImmutablePolynomial{M,R}, num)
    end
    deg = n - m +1
    q_coeff = zeros(R, deg+1)
    r_coeff = R[ num.a[i] for i in 1:n+1 ]
    @inbounds for i in n:-1:m
        q = r_coeff[i + 1] / den[m]
        q_coeff[i - m + 1] = q
        @inbounds for j in 0:m
            elem = den[j] * q
            r_coeff[i - m + j + 1] -= elem
        end
    end
    return ImmutablePolynomial(q_coeff, num.var), ImmutablePolynomial(r_coeff, num.var)
end
