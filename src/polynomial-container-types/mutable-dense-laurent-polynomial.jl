"""
    MutableDenseLaurentPolynomial{B,T,X}

This polynomial type essentially uses an offset vector (`Vector{T}`,`order`) to store the coefficients of a polynomial relative to the basis `B` with indeterminate `X`.

The typical offset is to have `0` as the order, but, say, to accomodate Laurent polynomials, or more efficient storage of basis elements any order may be specified.

This type trims trailing zeros and the leading zeros on construction.

"""
struct MutableDenseLaurentPolynomial{B,T,X} <: AbstractLaurentUnivariatePolynomial{B,T, X}
    coeffs::Vector{T}
    order::Base.RefValue{Int} # lowest degree, typically 0
    function MutableDenseLaurentPolynomial{B,T,X}(::Val{:false}, cs::AbstractVector, order::Int=0) where {B,T,X}
        new{B,T,Symbol(X)}(cs, Ref(order))
    end
    function MutableDenseLaurentPolynomial{B,T,X}(::Val{true}, cs::AbstractVector, order::Int=0) where {B,T,X}
        if Base.has_offset_axes(cs)
            @warn "Using the axis offset of the coefficient vector"
            cs, order = cs.parent, firstindex(cs)
        end

        i = findlast(!iszero, cs)
        if i == nothing
            xs = T[]
        else
            j = findfirst(!iszero, cs)
            xs = T[cs[i] for i ∈ j:i]
            order = order + j - 1
        end
        new{B,T,Symbol(X)}(xs, Ref(order))
    end
end

function MutableDenseLaurentPolynomial{B,T,X}(cs::AbstractVector{T}, order::Int=0) where {B,T,X}
    MutableDenseLaurentPolynomial{B,T,X}(Val(true), cs, order)
end

function _polynomial(p::P, as::Vector{S})  where {B,T, X, P <: MutableDenseLaurentPolynomial{B,T,X}, S}
    R = eltype(as)
    Q = MutableDenseLaurentPolynomial{B, R, X}
    as = trim_trailing_zeros!!(as)
    Q(Val(false), as, p.order[])
end

@poly_register MutableDenseLaurentPolynomial

## ---

## Generics for polynomials
function Base.convert(::Type{MutableDenseLaurentPolynomial{B,T,X}}, q::MutableDenseLaurentPolynomial{B,T′,X′}) where {B,T,T′,X,X′}
    MutableDenseLaurentPolynomial{B,T,X}(Val(false), convert(Vector{T},q.coeffs), q.order[])
end

function Base.map(fn, p::P, args...)  where {B,T,X, P<:MutableDenseLaurentPolynomial{B,T,X}}
    xs = map(fn, p.coeffs, args...)
    xs = trim_trailing_zeros!!(xs)
    R = eltype(xs)
    return ⟒(P){R,X}(Val(false), xs, p.order[])
end

Base.copy(p::MutableDenseLaurentPolynomial{B,T,X}) where {B,T,X} =
    MutableDenseLaurentPolynomial{B,T,Var(X)}(copy(p.coeffs), firstindex(p))

Base.firstindex(p::MutableDenseLaurentPolynomial) = p.order[]
Base.length(p::MutableDenseLaurentPolynomial) = length(p.coeffs)
Base.lastindex(p::MutableDenseLaurentPolynomial) = firstindex(p) + length(p) - 1
function Base.getindex(p::MutableDenseLaurentPolynomial{B,T,X}, i::Int) where {B,T,X}
    (i < firstindex(p) || i > lastindex(p)) && return zero(T)
    p.coeffs[i + offset(p)]
end
# ??? should this call chop! if `iszero(value)`?
function Base.setindex!(p::P, value, i::Int) where {B,T,X,P<:MutableDenseLaurentPolynomial{B,T,X}}
    a,b = firstindex(p), lastindex(p)
    o = a
    iszero(p) && return P([value], i)
    z = zero(first(p.coeffs))
    if i > b
        append!(p.coeffs, _zeros(p, z, i - b))
        p.coeffs[end] = value
    elseif i < a
        prepend!(p.coeffs, _zeros(p, z, a-i))
        p.coeffs[1] = value
        o = i
    else
        p.coeffs[i + offset(p)] = value
    end
    p.order[] = o
    p
end

offset(p::MutableDenseLaurentPolynomial) = 1 - firstindex(p)
Base.eachindex(p::MutableDenseLaurentPolynomial) = firstindex(p):1:lastindex(p)
Base.iterate(p::MutableDenseLaurentPolynomial, args...) = Base.iterate(p.coeffs, args...)
Base.pairs(p::MutableDenseLaurentPolynomial) =
    Base.Generator(=>, firstindex(p):lastindex(p), p.coeffs)

# return a container of zeros based on basis type
_zeros(::Type{<:MutableDenseLaurentPolynomial}, z, N)  = fill(z, N)

Base.similar(p::MutableDenseLaurentPolynomial, args...) = similar(p.coeffs, args...)

# iszero
Base.iszero(p::MutableDenseLaurentPolynomial) = iszero(p.coeffs)::Bool

function degree(p::MutableDenseLaurentPolynomial)
    i = findlast(!iszero, p.coeffs)
    isnothing(i) && return -1
    firstindex(p) + i - 1
end



basis(::Type{MutableDenseLaurentPolynomial{B,T,X}},i::Int) where {B,T,X} = MutableDenseLaurentPolynomial{B,T,X}([1],i)


function chop_left_index(x; rtol=nothing, atol=nothing)
    isempty(x) && return nothing
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, norm(x,2) * δ)
    i = findfirst(Base.Fix2(gtτ,τ), x)
    i
end

Base.chop(p::MutableDenseLaurentPolynomial{B,T,X}; kwargs...) where {B,T,X} = chop!(copy(p);kwargs...)

# trims left **and right**
function chop!(p::MutableDenseLaurentPolynomial{B,T,X};
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

function normΔ(q1::MutableDenseLaurentPolynomial{B}, q2::MutableDenseLaurentPolynomial{B}) where {B}
    iszero(q1) && return norm(q2,2)
    iszero(q2) && return norm(q1,2)
    r = abs(zero(q1[end] + q2[end]))
    tot = zero(r)
    for i ∈ minimum(firstindex,(q1,q2)):maximum(lastindex, (q1,q2))
       @inbounds tot += abs2(q1[i] - q2[i])
    end
    return sqrt(tot)
end

minimumexponent(::Type{<:MutableDenseLaurentPolynomial}) =  typemin(Int)



# vector ops +, -, c*x
## unary - (map is as fast)
## binary +
Base.:+(p::MutableDenseLaurentPolynomial{B,T,X}, q::MutableDenseLaurentPolynomial{B,S,X}) where{B,X,T,S} =
    offset_vector_combine(+, p, q)
Base.:-(p::MutableDenseLaurentPolynomial{B,T,X}, q::MutableDenseLaurentPolynomial{B,S,X}) where{B,X,T,S} =
    offset_vector_combine(-, p, q)
# modified from  https://github.com/jmichel7/LaurentPolynomials.jl/
function offset_vector_combine(op, p::MutableDenseLaurentPolynomial{B,T,X}, q::MutableDenseLaurentPolynomial{B,S,X}) where{B,X,T,S}
    R = promote_type(T,S)
    P = MutableDenseLaurentPolynomial{B,R,X}

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

    b₁ == b₂ && (x = trim_trailing_zeros!!(x))
    P(Val(false), x, a)

end

function Base.numerator(q::MutableDenseLaurentPolynomial{B,T,X}) where {B,T,X}
    p = chop(q)
    o = firstindex(p)
    o ≥ 0 && return p
    MutableDenseLaurentPolynomial{B,T,X}(p.coeffs, 0)
end

function Base.denominator(q::MutableDenseLaurentPolynomial{B,T,X}) where {B,T,X}
    p = chop(q)
    o = firstindex(p)
    o ≥ 0 && return one(p)
    basis(MutableDenseLaurentPolynomial{B,T,X}, -o)
end



## ---
function LinearAlgebra.lmul!(c::Scalar, p::MutableDenseLaurentPolynomial{B,T,X}) where {B,T,X}
    if iszero(c)
        empty!(p.coeffs)
        p.order[] = 0
    else
        lmul!(c, p.coeffs)
    end
    p
end

function LinearAlgebra.rmul!(p::MutableDenseLaurentPolynomial{B,T,X}, c::Scalar) where {B,T,X}
    if iszero(c)
        empty!(p.coeffs)
        p.order[] = 0
    else
        rmul!(p.coeffs, c)
    end
    p
end
