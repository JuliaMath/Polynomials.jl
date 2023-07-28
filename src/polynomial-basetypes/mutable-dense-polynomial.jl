"""
    MutableDensePolynomial{B,T,X}

This polynomial type essentially uses an offset vector (`Vector{T}`,`order`) to store the coefficients of a polynomial relative to the basis `B` with indeterminate `X`.

The typical offset is to have `0` as the order, but, say, to accomodate Laurent polynomials, or more efficient storage of basis elements any order may be specified.

This type trims trailing zeros and when the offset is not 0, trims the leading zeros.

"""
struct MutableDensePolynomial{B,T,X} <: AbstractUnivariatePolynomial{B,T, X}
    coeffs::Vector{T}
    order::Base.RefValue{Int} # lowest degree, typically 0
    function MutableDensePolynomial{B,T,X}(cs::AbstractVector{S}, order::Int=0) where {B,T,X,S}
        if Base.has_offset_axes(cs)
            @warn "Using the axis offset of the coefficient vector"
            cs, order = cs.parent, firstindex(cs)
        end

        i = findlast(!iszero, cs)
        if i == nothing
            xs = T[]
        elseif iszero(order)
            xs = T[cs[i] for i ∈ 1:i]
        else # shift if not 0
            j = findfirst(!iszero, cs)
            xs = T[cs[i] for i ∈ j:i]
            order = order + j - 1
        end
        new{B,T,Symbol(X)}(xs, Ref(order))

    end
    function MutableDensePolynomial{B,T,X}(check::Val{false}, cs::Vector{T}, order::Int=0) where {B,T,X}
        if Base.has_offset_axes(cs)
            @warn "Using the axis offset of the coefficient vector"
            cs, order = cs.parent, first(cs.offsets)
        end
        new{B,T,Symbol(X)}(cs, Ref(order))
    end
    function MutableDensePolynomial{B,T,X}(check::Val{true}, cs::Vector{T}, order::Int=0) where {B,T,X}
        MutableDensePolynomial{B,T,X}(cs, Ref(order))
    end
end

function _polynomial(p::P, as::Vector{S})  where {B,T, X, P <: MutableDensePolynomial{B,T,X}, S}
    R = eltype(as)
    Q = MutableDensePolynomial{B, R, X}
    as = trim_trailing_zeros(as)
    Q(Val(false), as, p.order[])
end

@poly_register MutableDensePolynomial
constructorof(::Type{<:MutableDensePolynomial{B}}) where {B} = MutableDensePolynomial{B}

## ---

## Generics for polynomials
function Base.convert(::Type{MutableDensePolynomial{B,T,X}}, q::MutableDensePolynomial{B,T′,X′}) where {B,T,T′,X,X′}
    MutableDensePolynomial{B,T,X}(Val(false), convert(Vector{T},q.coeffs), q.order[])
end

Base.copy(p::MutableDensePolynomial{B,T,X}) where {B,T,X} =
    MutableDensePolynomial{B,T,Var(X)}(copy(p.coeffs), firstindex(p))

Base.firstindex(p::MutableDensePolynomial) = p.order[]
Base.length(p::MutableDensePolynomial) = length(p.coeffs)
Base.lastindex(p::MutableDensePolynomial) = firstindex(p) + length(p) - 1
function Base.getindex(p::MutableDensePolynomial{B,T,X}, i::Int) where {B,T,X}
    (i < firstindex(p) || i > lastindex(p)) && return zero(T)
    p.coeffs[i + offset(p)]
end
# ??? should this call chop! if `iszero(value)`?
function Base.setindex!(p::P, value, i::Int) where {B,T,X,P<:MutableDensePolynomial{B,T,X}}
    a,b = firstindex(p), lastindex(p)
    o = a
    iszero(p) && return P([value], i)
    z = zero(first(p.coeffs))
    if i > b
        append!(p.coeffs, _zeros(p, z, i - b))
        p.coeffs[end] = value
    elseif i < a
        prepend!(p.coeffs, zeros(p, z, a-i))
        p.coeffs[1] = value
        o = i
    else
        p.coeffs[i + offset(p)] = value
    end
    p.order[] = o
    p
end

offset(p::MutableDensePolynomial) = 1 - firstindex(p)
Base.eachindex(p::MutableDensePolynomial) = firstindex(p):1:lastindex(p)
Base.iterate(p::MutableDensePolynomial, args...) = Base.iterate(p.coeffs, args...)
Base.pairs(p::MutableDensePolynomial) =
    Base.Generator(=>, firstindex(p):lastindex(p), p.coeffs)

# return a container of zeros based on basis type
_zeros(::Type{<:MutableDensePolynomial}, z, N)  = fill(z, N)

Base.similar(p::MutableDensePolynomial, args...) = similar(p.coeffs, args...)

# iszero
Base.iszero(p::MutableDensePolynomial) = iszero(p.coeffs)::Bool

function degree(p::MutableDensePolynomial)
    i = findlast(!iszero, p.coeffs)
    isnothing(i) && return -1
    firstindex(p) + i - 1
end


laurenttype(::Type{<:MutableDensePolynomial}) = Val(true)

basis(::Type{MutableDensePolynomial{B,T,X}},i::Int) where {B,T,X} = MutableDensePolynomial{B,T,X}([1],i)

function trim_trailing_zeros(cs::Vector{T}) where {T}
    isempty(cs) && return cs
    !iszero(last(cs)) && return cs
    i = findlast(!iszero, cs)
    if isnothing(i)
        empty!(cs)
    else
        n = length(cs)
        while n > i
            pop!(cs)
            n -= 1
        end
    end
    cs
end



Base.chop(p::MutableDensePolynomial{B,T,X}; kwargs...) where {B,T,X} = chop!(copy(p);kwargs...)

# trims left **and right**
function chop!(p::MutableDensePolynomial{B,T,X};
               atol=nothing, rtol=Base.rtoldefault(float(real(T)))) where {B,T,X}
    iᵣ = chop_right_index(p.coeffs; atol=atol, rtol=rtol) # nothing -- nothing to chop
    iᵣ === nothing && return zero(p)
    iₗ = chop_left_index(p.coeffs; atol=atol, rtol=rtol)
    iₗ === nothing && return zero(p)
    iₗ == 1 && iᵣ == length(p.coeffs) && return p

    N = length(p.coeffs)

    o = firstindex(p)
    for i ∈ 1:(iₗ-1)
        popfirst!(p.coeffs)
        o += 1
    end

    for i ∈ (iᵣ+1):N
        pop!(p.coeffs)
    end
    p.order[] = o
    p
end

function normΔ(q1::MutableDensePolynomial{B}, q2::MutableDensePolynomial{B}, p::Real = 2) where {B}
    iszero(q1) && return norm(q2, p)
    iszero(q2) && return norm(q1, p)
    r = zero(q1[end] + q2[end])
    tot = zero(r)
    for i ∈ minimum(firstindex,(q1,q2)):maximum(lastindex, (q1,q2))
       @inbounds tot += (q1[i] - q2[i])^p
    end
    return tot^(1/p)
end

minimumexponent(::Type{<:MutableDensePolynomial}) =  typemin(Int)



# vector ops +, -, c*x
## unary
Base.:-(p::MutableDensePolynomial{B,T,X}) where {B,T,X} =
    MutableDensePolynomial{B,T,X}(Val(false), -p.coeffs, p.order[])

## binary
Base.:+(p::MutableDensePolynomial{B,T,X}, q::MutableDensePolynomial{B,S,X}) where{B,X,T,S} =
    offset_vector_combine(+, p, q)
Base.:-(p::MutableDensePolynomial{B,T,X}, q::MutableDensePolynomial{B,S,X}) where{B,X,T,S} =
    offset_vector_combine(-, p, q)
# handle +, -
# modified from  https://github.com/jmichel7/LaurentPolynomials.jl/
function offset_vector_combine(op, p::MutableDensePolynomial{B,T,X}, q::MutableDensePolynomial{B,S,X}) where{B,X,T,S}
    R = promote_type(T,S)
    P = MutableDensePolynomial{B,R,X}

    iszero(p) && return convert(P, op(q))
    iszero(q) && return convert(P, p)

    a₁, a₂ = firstindex(p), firstindex(q)
    b₁, b₂ = lastindex(p), lastindex(q)
    a, b = min(a₁, a₂), max(b₁, b₂)

    N = b - a + 1
    z = zero(first(p.coeffs) + first(q.coeffs))
    x = _zeros(p, z, N)

    Δp, Δq =  a₁ - a₂, 0
    if a₁ < a₂
        Δq, Δp = -Δp, Δq
    end
    # zip faster than `pairs`
    @inbounds for (i, cᵢ) ∈ zip((1+Δp):(length(p) + Δp), p.coeffs)
        x[i] = cᵢ
    end
    @inbounds for (i, cᵢ) ∈ zip((1+Δq):(length(q) + Δq), q.coeffs)
        x[i] = op(x[i], cᵢ)
    end

    b₁ == b₂ && (x = trim_trailing_zeros(x))
    P(Val(false), x, a)

end


function Base.numerator(p::MutableDensePolynomial{B,T,X}) where {B,T,X}
    o = firstindex(p)
    o ≥ 0 && return p
    MutableDensePolynomial{B,T,X}(p.coeffs, 0)
end

function Base.denominator(p::MutableDensePolynomial{B,T,X}) where {B,T,X}
    o = firstindex(p)
    o ≥ 0 && return one(p)
    basis(MutableDensePolynomial{B,T,X}, -o)
end


# scalar mult
function scalar_mult(p::MutableDensePolynomial{B,T,X}, c::S) where {B,T,X,S}
    cs = p.coeffs .* (c,) # works with T[]
    return _polynomial(p, cs)
end
function scalar_mult(c::S, p::MutableDensePolynomial{B,T,X}) where {B,T,X,S}
    cs = (c,) .* p.coeffs
    return _polynomial(p, cs)
end

## ---
function LinearAlgebra.lmul!(c::Scalar, p::MutableDensePolynomial{B,T,X}) where {B,T,X}
    if iszero(c)
        empty!(p.coeffs)
        p.order[] = 0
    else
        lmul!(c, p.coeffs)
    end
    p
end

function LinearAlgebra.rmul!(p::MutableDensePolynomial{B,T,X}, c::Scalar) where {B,T,X}
    if iszero(c)
        empty!(p.coeffs)
        p.order[] = 0
    else
        rmul!(p.coeffs, c)
    end
    p
end
