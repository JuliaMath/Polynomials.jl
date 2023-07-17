struct StandardBasis <: AbstractBasis end

function showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {B<:StandardBasis, T, P<:AbstractUnivariatePolynomial{B,T}}

    if _iszero(pj) return false end

    pj = printsign(io, pj, first, mimetype)

    if hasone(T)
        if !(_isone(pj) && !(showone(T) || j == 0))
            printcoefficient(io, pj, j, mimetype)
        end
    else
        printcoefficient(io, pj, j, mimetype)
    end

    printproductsign(io, pj, j, mimetype)
    printexponent(io, var, j, mimetype)
    return true
end


function print_basis(io::IO, p::AbstractUnivariatePolynomial{<:StandardBasis, T, X}, i) where {T,X}
    print(io, X)
    print_unicode_exponent(io, i)
end

function Base.one(::Type{P}) where {B<:StandardBasis,T,X, P <: AbstractUnivariatePolynomial{B,T,X}}
    ⟒(P){B,T,X}([1])
end

function variable(P::Type{<:AbstractUnivariatePolynomial{B,T,X}}) where {B<:StandardBasis,T,X}
    basis(P, 1)
end

function basis(P::Type{<:AbstractUnivariatePolynomial{B, T, X}}, i) where {B<:StandardBasis,T,X}
    cs = ones(T,1)
    P(cs, i)
end

constantterm(p::AbstractUnivariatePolynomial{B}) where {B <: StandardBasis} = p[0]



# # storage independent scalar add
# # faster alternatives for basic types
# function scalar_add(c::S, p::AbstractUnivariatePolynomial{B,T,X}) where {B<:StandardBasis, S, T, X}
#     R = promote_type(T,S)
#     P = ⟒(p){B,R,X}

#     iszero(p) && return P((c,), 0)
#     iszero(c) && return convert(P, p)

#     a,b = firstindex(p), lastindex(p)
#     a′ = min(0,a)
#     z = zero(last(first(p.coeffs)) + c)
#     cs = _zeros(p, zero(z), length(a′:b))

#     # i -> idx mapping
#     o = offset(p)*(-a+1) # offset for dict is 0; o/w -a + 1
#     for (i, cᵢ) ∈ pairs(p)
#         cs = _set(cs, i+o, R(cᵢ))
#     end
#     cs = _set(cs, 0+o, cs[0+o] + R(c))
#     cs = trim_trailing_zeros(cs)
#     P(Val(false), cs, a′)
# end

# special cases are faster
function ⊗(p::AbstractUnivariatePolynomial{B,T,X},
           q::AbstractUnivariatePolynomial{B,S,X}) where {B <: StandardBasis, T,S,X}
    # simple convolution with order shifted
    R = promote_type(T,S)
    P = ⟒(p){B,R,X}

    iszero(p) && return zero(P)
    iszero(q) && return zero(P)

    cs = fastconv(p.coeffs, q.coeffs)
    R = eltype(P)
    a = firstindex(p) + firstindex(q)

    P(cs, a)
end

# maybe do same for, say, Laguerre? Though that may push everything to Dense
function differentiate(p::AbstractUnivariatePolynomial{B,T,X}) where {B<:StandardBasis,T,X}
    N = lastindex(p) - firstindex(p) + 1
    R = promote_type(T, Int)
    P = ⟒(p){B,T,X}
    iszero(p) && return zero(P)
    z = zero(1 * p[1])
    cs = _zeros(p, z, N)
    os = offset(p)
    @inbounds for (i, cᵢ) ∈ pairs(p)
        iszero(i) && continue
        #cs[i - 1 + os] = i * cᵢ
        cs = _set(cs, i - 1 + os, i * cᵢ)
    end

    o = firstindex(p)
    o = o < 0 ? o - 1 : max(0, o - 1)
    ⟒(p){B,T,X}(cs, o)
end



function integrate(p::AbstractUnivariatePolynomial{B,T,X}) where {B <: StandardBasis,T,X}
    # no offset! XXX

    iszero(p) && return p/1

    N = lastindex(p) - firstindex(p) + 1
    R = typeof(one(T)/1)
    z = zero(R)
    P = ⟒(p){B,R,X}
    cs = _zeros(p, z, N+1)
    os =  offset(p)
    @inbounds for (i, cᵢ) ∈ pairs(p)
        i == -1 && throw(ArgumentError("Laurent polynomial with 1/x term"))
        #cs[i + os] = cᵢ / (i+1)
        cs = _set(cs, i + 1 + os,  cᵢ / (i+1))
    end
    P(cs, firstindex(p))
end
