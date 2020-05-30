abstract type StandardBasisPolynomial{T} <: AbstractPolynomial{T} end



function showterm(io::IO, ::Type{<:StandardBasisPolynomial}, pj::T, var, j, first::Bool, mimetype) where {T} 
    if iszero(pj) return false end
    pj = printsign(io, pj, first, mimetype)
    if !(isone(pj) && !(showone(T) || iszero(j)))
        printcoefficient(io, pj, j, mimetype)
    end
    printproductsign(io, pj, j, mimetype)
    printexponent(io, var, j, mimetype)
    return true
end

# allows  broadcast  issue #209
evalpoly(x, p::StandardBasisPolynomial) = p(x)

domain(::Type{<:StandardBasisPolynomial}) = Interval(-Inf, Inf)
mapdomain(::Type{<:StandardBasisPolynomial}, x::AbstractArray) = x

## generic test if polynomial `p` is a constant
isconstant(p::StandardBasisPolynomial) = degree(p) <= 0

Base.convert(P::Type{<:StandardBasisPolynomial}, q::StandardBasisPolynomial) = isa(q, P) ? q : P([q[i] for i in 0:degree(q)], q.var)

variable(::Type{P}, var::SymbolLike = :x) where {P <: StandardBasisPolynomial} = P([0, 1], var)

function fromroots(P::Type{<:StandardBasisPolynomial}, r::AbstractVector{T}; var::SymbolLike = :x) where {T <: Number}
    n = length(r)
    c = zeros(T, n + 1)
    c[1] = one(T)
    for j in 1:n, i in j:-1:1
        c[(i + 1)] = c[(i + 1)] - r[j] * c[i]
    end
    #return P(c, var)
    return P(reverse(c), var)
end


function Base.:+(p::P, c::S) where {T, P <: StandardBasisPolynomial{T}, S}
    R = promote_type(T,S)
    as = R[c  for c in coeffs(p)]
    as[1] += c
    ⟒(P)(as, p.var)
end


function derivative(p::P, order::Integer = 1) where {T, P <: StandardBasisPolynomial{T}}
    order < 0 && error("Order of derivative must be non-negative")

    # we avoid usage like Base.promote_op(*, T, Int) here, say, as
    # Base.promote_op(*, Rational, Int) is Any, not Rational in analogy to
    # Base.promote_op(*, Complex, Int)
    R = eltype(one(T)*1) 
    order == 0 && return p
    hasnan(p) && return ⟒(P){R}(R[NaN], p.var)
    order > length(p) && return zero(⟒(P){R},p.var)
    d = degree(p)
    d <= 0 && return zero(⟒(P){R},p.var)
    n = d + 1
    a2 = Vector{R}(undef, n - order)
    @inbounds for i in order:n - 1
        a2[i - order + 1] = reduce(*, (i - order + 1):i, init = p[i])
    end
    return ⟒(P)(a2, p.var)
end


function integrate(p::P, k::S) where {T, P <: StandardBasisPolynomial{T}, S<:Number}

    R = eltype((one(T)+one(S))/1)
    if hasnan(p) || isnan(k)
        return ⟒(P)([NaN])
    end
    n = length(p)
    a2 = Vector{R}(undef, n + 1)
    a2[1] = k
    @inbounds for i in 1:n
        a2[i + 1] = p[i - 1] / i
    end
    return ⟒(P)(a2, p.var)
end


function Base.divrem(num::P, den::Q) where {T, P <: StandardBasisPolynomial{T}, S, Q <: StandardBasisPolynomial{S}}

    check_same_variable(num, den) || error("Polynomials must have same variable")
    var = num.var
    
    
    n = degree(num)
    m = degree(den)

    m == -1 && throw(DivideError())
    if m == 0 && den[0] ≈ 0 throw(DivideError()) end
    
    R = eltype(one(T)/one(S))

    deg = n - m + 1

    if deg ≤ 0
        return zero(P, var), num
    end

    q_coeff = zeros(R, deg)
    r_coeff = R[ num[i-1] for i in 1:n+1 ]

    @inbounds for i in n:-1:m
        q = r_coeff[i + 1] / den[m]
        q_coeff[i - m + 1] = q
        @inbounds for j in 0:m
            elem = den[j] * q
            r_coeff[i - m + j + 1] -= elem
        end
    end

    return ⟒(P)(q_coeff, var), ⟒(P)(r_coeff, var)
    
end


function companion(p::P) where {T, P <: StandardBasisPolynomial{T}}
    d = length(p) - 1
    d < 1 && error("Series must have degree greater than 1")
    d == 1 && return diagm(0 => [-p[0] / p[1]])

    
    R = eltype(one(T)/one(T))
    
    comp = diagm(-1 => ones(R, d - 1))
    ani = 1 / p[end]
    for j in  0:(degree(p)-1)
        comp[1,(d-j)] = -p[j] * ani # along top row has smaller residual than down column
    end
    return comp
end

function  roots(p::P; kwargs...)  where  {T, P <: StandardBasisPolynomial{T}}

    R = eltype(one(T)/one(T))    
    d = degree(p)
    if d < 1
        return []
    end
    d == 1 && return R[-p[0] / p[1]]

    as = [p[i] for i in 0:d]
    K  = findlast(!iszero, as)
    if K == nothing
        return R[]
    end
    k =  findfirst(!iszero, as)

    k  == K && return zeros(R, k-1)

    # find eigenvalues of the  companion matrix
    comp  = companion(⟒(P)(as[k:K], p.var))
    L = eigvals(comp; kwargs...)
    append!(L, zeros(eltype(L), k-1))

    L
end


function vander(P::Type{<:StandardBasisPolynomial}, x::AbstractVector{T}, n::Integer) where {T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= one(T)
    @inbounds for i in 1:n
        A[:, i + 1] = A[:, i] .* x
    end
    return A
end

