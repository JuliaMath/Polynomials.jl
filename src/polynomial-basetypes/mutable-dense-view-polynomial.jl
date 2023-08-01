"""
    MutableDenseViewPolynomial{B,T,X}

Construct a polynomial in `P_n` (or `Πₙ`), the collection of polynomials in the
basis of degree `n` *or less*, using a vector of length
`N+1`.

* Unlike other polynomial types, this type allows trailing zeros in the coefficient vector. Call `chop!` to trim trailing zeros if desired.
* Unlike other polynomial types, this does not copy the coefficients on construction
* Unlike other polynomial types, this type broadcasts like a vector for in-place vector operations (scalar multiplication, polynomial addition/subtraction of the same size)
* The method inplace `mul!(pq, p, q)` is defined to use precallocated storage for the product of `p` and `q`

This type is useful for reducing copies and allocations in some algorithms.

"""
struct MutableDenseViewPolynomial{B,T,X} <: AbstractUnivariatePolynomial{B, T, X}
    coeffs::Vector{T}
    function MutableDenseViewPolynomial{B, T, X}(coeffs::AbstractVector{S}) where {B, T,S, X}
        new{B,T,Symbol(X)}(convert(Vector{T}, coeffs))
    end
end

MutableDensePolynomial{B,T,X}(check::Val{false}, cs::AbstractVector{S}) where {B,T,S,X} = new{B,T,X}(coeffs)
MutableDensePolynomial{B,T,X}(check::Val{true}, cs::AbstractVector{S}) where {B,T,S,X} = new{B,T,X}(coeffs)

function MutableDenseViewPolynomial{B, T, X}(coeffs::AbstractVector{S}, order::Int) where {B, T,S, X}
    iszero(order) && return MutableDenseViewPolynomial{B,T,X}(coeffs)
    order < 0 && throw(ArgumentError("Not a Laurent type"))
    prepend!(coeffs, zeros(T,order))
    MutableDenseViewPolynomial{B,T,X}(coeffs)
end


@poly_register MutableDenseViewPolynomial
constructorof(::Type{<:MutableDenseViewPolynomial{B}}) where {B} = MutableDenseViewPolynomial{B}
minimumexponent(::Type{<:MutableDenseViewPolynomial}) =  0
laurenttype(::Type{<:MutableDenseViewPolynomial}) = Val(false)
_zeros(::Type{<:MutableDenseViewPolynomial}, z, N)  = fill(z, N)
Base.copy(p::MutableDenseViewPolynomial{B,T,X}) where {B,T,X} = MutableDenseViewPolynomial{B,T,X}(copy(p.coeffs))
# change broadcast semantics
Base.broadcastable(p::MutableDenseViewPolynomial) = p.coeffs;
Base.ndims(::Type{<:MutableDenseViewPolynomial}) = 1
Base.copyto!(p::MutableDenseViewPolynomial{B, T, X},
             x::S) where {B, T, X,
                          S<:Union{AbstractVector, Base.AbstractBroadcasted, Tuple} # to avoid an invalidation. Might need to be more general?
                          } = copyto!(p.coeffs, x)

Base.firstindex(p::MutableDenseViewPolynomial) = 0
Base.lastindex(p::MutableDenseViewPolynomial) = length(p.coeffs)-1
Base.pairs(p::MutableDenseViewPolynomial) = Base.Generator(=>, 0:(length(p.coeffs)-1), p.coeffs)
function degree(p::MutableDenseViewPolynomial)
    i = findlast(!iszero, p.coeffs)
    isnothing(i) && return -1
    i - 1
end

# trims left **and right**
function chop!(p::MutableDenseViewPolynomial{B,T,X};
               atol=nothing, rtol=Base.rtoldefault(float(real(T)))) where {B,T,X}
    iᵣ = chop_right_index(p.coeffs; atol=atol, rtol=rtol) # nothing -- nothing to chop
    iᵣ === nothing && return zero(p)
    N = length(p.coeffs)
    for i ∈ (iᵣ+1):N
        pop!(p.coeffs)
    end
    p
end

function Base.:-(p::MutableDenseViewPolynomial{B,T,X}) where {B,T,X}
    MutableDenseViewPolynomial{B,T,X}(-p.coeffs) # use lmul!(-1,p) for in place
end

# for same length, can use p .+= q for in place
Base.:+(p::MutableDenseViewPolynomial{B,T,X}, q::MutableDenseViewPolynomial{B,S,X}) where{B,X,T,S} =
    _vector_combine(+, p, q)
Base.:-(p::MutableDenseViewPolynomial{B,T,X}, q::MutableDenseViewPolynomial{B,S,X}) where{B,X,T,S} =
    _vector_combine(-, p, q)

function _vector_combine(op, p::MutableDenseViewPolynomial{B,T,X}, q::MutableDenseViewPolynomial{B,S,X}) where {B,T,S,X}
    n,m = length(p.coeffs), length(q.coeffs)
    R = promote_type(T,S)
    if n ≥ m
        cs = convert(Vector{R}, p.coeffs)
        for (i,qᵢ) ∈ enumerate(q.coeffs)
            cs[i] = op(cs[i], qᵢ)
        end
    else
        cs = convert(Vector{R}, q.coeffs)
        for (i,pᵢ) ∈ enumerate(p.coeffs)
            cs[i] = op(pᵢ, cs[i])
        end
    end
    MutableDenseViewPolynomial{B,R,X}(cs)
end


# pre-allocated multiplication
function LinearAlgebra.lmul!(c::Number, p::MutableDenseViewPolynomial{T,X}) where {T,X}
    p.coeffs[:] = (c,) .* p.coeffs
    p
end
function LinearAlgebra.rmul!(p::MutableDenseViewPolynomial{T,X}, c::Number) where {T,X}
    p.coeffs[:] = p.coeffs .* (c,)
    p
end
