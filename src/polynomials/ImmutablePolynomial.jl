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
    function ImmutablePolynomial{N, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {N, T <: Number}
        length(coeffs) == 0 && return new{1, T}(tuple(zeros(T, 1)...), var)
        last_nz = findlast(!iszero, coeffs)
        M = max(1, last_nz === nothing ? 0 : last_nz)
        M > N &&  throw(ArgumentError(""))
        return new{N,T}(NTuple{N, T}(i <= M ? coeffs[i] : zero(T) for i in  1:N), var)
    end
end

## Can't do this so quickly, as we keep the parameter N
## @register Polynomial
Base.convert(::Type{P}, p::P) where {P <: ImmutablePolynomial} = p
function Base.convert(::Type{ImmutablePolynomial{N,T}}, p::ImmutablePolynomial{M,S}) where {N,T,M,S}
    N >= M || throw(ArgumentError("XXX"))
    ImmutablePolynomial{N,T}(NTuple{N,T}(T(p[i]) for i in 0:N-1), p.var)
end
Base.convert(P::Type{<:ImmutablePolynomial}, p::ImmutablePolynomial)  = P(coeffs(p), p.var)
Base.convert(P::Type{<:ImmutablePolynomial}, p::Polynomial)  = P(coeffs(p), p.var)
Base.convert(P::Type{<:Polynomial}, p::ImmutablePolynomial{N,T}) where {N,T}  = P(T[coeffs(p)...], p.var)
Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{ImmutablePolynomial{N,S}}) where  {N, T,S} =
    ImmutablePolynomial{N, promote_type(T, S)}
Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{ImmutablePolynomial{M,S}}) where  {N, T,M,S} =
    ImmutablePolynomial{max(N,M), promote_type(T, S)}
Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{P}) where {N,T, P<:AbstractPolynomial} =
    P
Base.promote_rule(::Type{ImmutablePolynomial{N,T}}, ::Type{S}) where  {N,T,S<:Number}  =
    ImmutablePolynomial{N, promote_type(T, S)}
function ImmutablePolynomial{N,T}(coeffs::NTuple{M,S}, var::Polynomials.SymbolLike = :x) where{N,T, M,S} 
    N >= M  || throw(ArgumentError(""))
    ImmutablePolynomial{N,T}(NTuple{N, T}(i <= M ? T(coeffs[i]) : zero(T) for i in 1:N), Symbol(var))
end
ImmutablePolynomial{N}(coeffs::NTuple{M,T}, var::Polynomials.SymbolLike = :x) where{N,T, M} =
    ImmutablePolynomial{M,T}(coeffs, Symbol(var))
function ImmutablePolynomial{N}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where{N,T}
    M = length(coeffs)
    N >= M || throw(ArgumentError(""))
    ImmutablePolynomial{N,T}(NTuple{N,T}(i <= M ? T(coeffs[i]) : zero(T) for i in  1:N), Symbol(var))
end
ImmutablePolynomial(coeffs::NTuple{N,T}, var::Polynomials.SymbolLike = :x) where{N,T} =
    ImmutablePolynomial{N,T}(coeffs, Symbol(var))
function ImmutablePolynomial(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where{N,T}
    M = length(coeffs)
    ImmutablePolynomial{M,T}(NTuple{M,T}(x for x in coeffs), Symbol(var))
end
ImmutablePolynomial{N,T}(coeffs::AbstractVector{S}, var::Polynomials.SymbolLike = :x) where{N,T, S} =
    ImmutablePolynomial{N,T}(T.(coeffs), Symbol(var))
ImmutablePolynomial{N,T}(n::Number, var::Polynomials.SymbolLike = :x) where {N,T} =
    ImmutablePolynomial{N,T}(NTuple{N,T}(i==1 ? T(n) : zero(T) for i in 1:N), Symbol(var))
ImmutablePolynomial{N}(n::T, var::Polynomials.SymbolLike = :x) where {N, T <: Number} =
    ImmutablePolynomial{N,T}(NTuple{N,T}(i==1 ? n : zero(T) for i in 1:N), Symbol(var))
ImmutablePolynomial{N,T}(var::Polynomials.SymbolLike = :x) where {N, T} = variable(ImmutablePolynomial{N,T}, var)
ImmutablePolynomial{N}(var::Polynomials.SymbolLike = :x) where {N} = variable(ImmutablePolynomial{N}, var)
ImmutablePolynomial(var::Polynomials.SymbolLike = :x)  = variable(ImmutablePolynomial{1, Int}, var)
ImmutablePolynomial(n::S, var::Polynomials.SymbolLike = :x) where {S <: Number} =
    ImmutablePolynomial{1,S}(NTuple{1,S}(n))

# overrides from common.jl due to coeffs possibly being padded
Base.copy(p::P) where {P <: ImmutablePolynomial} = P(coeffs(p), p.var)
function Base.hash(p::ImmutablePolynomial{N,T}, h::UInt) where {N,T}
    n = findlast(!iszero, coeffs(p))
    isnothing(n) && return hash(p.var, hash(NTuple{0,T}(),h))
    hash(p.var, hash(coeffs(p)[1:n], h))
end
Base.isequal(p::ImmutablePolynomial, q::ImmutablePolynomial) = hash(p) == hash(q)
Base.one(::Type{ImmutablePolynomial{N,T}}) where {N,T}  = ImmutablePolynomial(tuple((i==1 ? one(T) : zero(T) for i in  1:N)...))
Base.one(::Type{<:ImmutablePolynomial{N}}) where {N}    = ImmutablePolynomial(tuple((i==1 ? 1 : 0 for i in  1:N)...))
Base.zero(::Type{ImmutablePolynomial{N,T}}) where {N,T} = ImmutablePolynomial(zeros(T, N))
Base.zero(::Type{<:ImmutablePolynomial{N}}) where {N}   = ImmutablePolynomial(zeros(Int,N))
function degree(p::ImmutablePolynomial{N,T}) where {N, T}
    n = findlast(!iszero, coeffs(p))
    isnothing(n) && return -1
    n-1
end
function Base.:(==)(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    (p1.var == p2.var) || return false
    p1s, p2s = coeffs(p1), coeffs(p2)
    (N == M  && p1s == p2s) &&  return  true
    n1 = findlast(!iszero, p1s)
    n2 = findlast(!iszero, p2s)
    isnothing(n1) && isnothing(n2) && return true
    (isnothing(n1) || isnothing(n2)) && return false
    p1s[1:n1] == p2s[1:n2] &&  return true
    false
end
function Base.:(==)(p::ImmutablePolynomial{N,T}, n::Number) where {N,T}
    p[0] == n  || return false
    for i in 1:N-1
        iszero(p[i]) || return false
    end
    true
end
    
function Base.isapprox(p1::ImmutablePolynomial{N,T},
                       p2::ImmutablePolynomial{M,S};
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {N,T,M,S}
    p1.var == p2.var || error("p1 and p2 must have same var")
    NN = max(N,M)
    for i in 1:NN-1
        isapprox(p1[i],p2[i], rtol=rtol, atol=atol) || return false
    end
    true
end
function Base.isapprox(p1::ImmutablePolynomial{N,T}, n::S;
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {N, T,S}
    N >= 1 || return false
    isapprox(p1[0], n, rtol=rtol, atol=atol) || return false
    for i in 1:N-1
        isapprox(p1[i],zero(T), rtol=rtol, atol=atol) || return false
    end
    return true
end
function Base.chop(p::ImmutablePolynomial{N,T};
              rtol::Real = Base.rtoldefault(real(T)),
              atol::Real = 0)  where {N, T}
    cs = coeffs(p)
    for i in N:-1:1
        if !isapprox(cs[i], zero(T), rtol=rtol, atol=atol)
            return ImmutablePolynomial{i,T}(cs[1:i], p.var)
        end
    end
    zero(ImmutablePolynomial{0,T})
end
function Base.truncate(p::ImmutablePolynomial{N,T};
                       rtol::Real = Base.rtoldefault(real(T)),
                       atol::Real = 0)  where {N, T}
    q = chop(p, rtol=rtol, atol=atol)
    iszero(q) && return  q
    cs = coeffs(q)
    thresh = maximum(abs,cs) * rtol + atol
    ImmutablePolynomial{length(cs), T}(map(c->abs(c) <= thresh ? zero(T) : c, coeffs(q)), p.var)
end    
LinearAlgebra.conj(p::P) where {P <: ImmutablePolynomial} = P(conj([coeffs(p)...]))

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

function integrate(p::ImmutablePolynomial{N, T}, c::Number) where {N, T}
    R = typeof(one(T)/1)
    q = integrate(convert(Polynomial, p), c)
    convert(ImmutablePolynomial{N+1,R}, q)
end

function derivative(p::ImmutablePolynomial{N,T}, order::Integer = 1) where {N,T}
    order < 0 && error("Order of derivative must be non-negative")
    R = eltype(p[0]/1)
    order == 0 && return convert(ImmutablePolynomial{N,R}, p)
    hasnan(p) && return ImmutablePolynomial(R[NaN], p.var)
    order > length(p) && return zero(ImmutablePolynomial{0,R})

    n = length(p)
    a2 = Vector{R}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return ImmutablePolynomial{n-order,R}(NTuple{n-order, R}(a   for  a in a2), p.var)
end

companion(p::ImmutablePolynomial{T}) where T = companion(convert(Polynomial{T},p))
roots(p::ImmutablePolynomial{N,T}; kwargs...)  where  {N,T} = roots(convert(Polynomial{T}, p); kwargs...)

function Base.:+(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    R = promote_type(T,S)
    m = max(N, M)
    return ImmutablePolynomial{m, R}(NTuple{m, R}(p1[i] + p2[i] for i in 0:m-1), p1.var)
end


function Base.:+(p::ImmutablePolynomial{N, T}, c::S) where {N, T,S<:Number}
    R = promote_type(T, S)
    cs = R[pi for pi in coeffs(p)]
    if isempty(cs)
        return ImmutablePolynomial{1,R}(NTuple{1,R}(R(c)), p.var)
    else
        cs[1] += c
    end
    ImmutablePolynomial{N,R}(NTuple{N,R}(c for c in cs), p.var)
end

function Base.:*(p1::ImmutablePolynomial{N,T}, p2::ImmutablePolynomial{M,S}) where {N,T,M,S}
    p1.var != p2.var && error("Polynomials must have same variable")
    N == 0  &&  return  p1; M ==  0  &&  return  p2
    R = promote_type(T, S)
    cs = zeros(R, N + M -1)
    for i in 0:N-1, j in 0:M-1
        @inbounds cs[i + j + 1] += p1[i] * p2[j]
    end
    return ImmutablePolynomial{N+M-1,R}(NTuple{N+M-1,R}(c for c in cs), p1.var)
end

function Base.:*(p::ImmutablePolynomial{N,T}, c::S) where {N,T,S <: Number}
    R = eltype(one(T)*one(S))
    return ImmutablePolynomial{N,R}(NTuple{N,R}(p[i]*c for i in eachindex(p)), p.var)
end

Base.:-(p::ImmutablePolynomial{N,T}) where {N,T} = ImmutablePolynomial(NTuple{N,T}(-pi for pi in p.coeffs), p.var)
Base.:*(c::Number, p::ImmutablePolynomial) = *(p, c)
Base.:/(p::ImmutablePolynomial, c::S) where {S} = ImmutablePolynomial([a/c for a in coeffs(p)])
Base.to_power_type(p::ImmutablePolynomial{N,T}) where {N,T} = p


function Base.divrem(num::ImmutablePolynomial{N, T}, den::ImmutablePolynomial{M,S}) where {N,T,M,S}
    num.var != den.var && error("Polynomials must have same variable")
    n = degree(num)
    m = degree(den)
    m == -1 && throw(DivideError())
    R = typeof(one(T) / one(S))
    if m > n
        return zero(ImmutablePolynomial{1,R}), ImmutablePolynomial{M,R}(num.coeffs, num.var)
    end
    deg = n - m +1
    q_coeff = zeros(R, deg+1)
    r_coeff = R[ coeffs(num)[i] for i in 1:n+1 ]
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
