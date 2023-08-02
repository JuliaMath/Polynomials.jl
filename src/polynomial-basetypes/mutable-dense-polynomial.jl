"""
    MutableDensePolynomial{B,T,X}

"""
struct MutableDensePolynomial{B,T,X} <: AbstractDenseUnivariatePolynomial{B,T, X}
    coeffs::Vector{T}
    function MutableDensePolynomial{B,T,X}(::Val{true}, cs::AbstractVector{S}, order::Int=0) where {B,T,X,S}
        if Base.has_offset_axes(cs)
            @warn "ignoring the axis offset of the coefficient vector"
            cs = parent(cs)
        end
        i = findlast(!iszero, cs)
        if i == nothing
            xs = T[]
        else
            xs = T[cs[i] for i ∈ 1:i] # make copy
        end
        if order > 0
            prepend!(xs, zeros(T, order))
        end
        new{B,T,Symbol(X)}(xs)

    end
    function MutableDensePolynomial{B,T,X}(check::Val{false}, cs::AbstractVector{T}, order::Int=0) where {B,T,X}
        new{B,T,Symbol(X)}(cs)
    end
end
function MutableDensePolynomial{B,T,X}(cs::AbstractVector{T}, order::Int=0) where {B,T,X}
    MutableDensePolynomial{B,T,X}(Val(true), cs, order)
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
    MutableDensePolynomial{B,T,X}(Val(false), convert(Vector{T},q.coeffs))
end

function Base.map(fn, p::P, args...)  where {B,T,X, P<:MutableDensePolynomial{B,T,X}}
    xs = map(fn, p.coeffs, args...)
    xs = trim_trailing_zeros(xs)
    R = eltype(xs)
    return ⟒(P){R,X}(Val(false), xs)
end

Base.copy(p::MutableDensePolynomial{B,T,X}) where {B,T,X} =
    MutableDensePolynomial{B,T,Var(X)}(copy(p.coeffs))

Base.firstindex(p::MutableDensePolynomial) = 0
Base.length(p::MutableDensePolynomial) = length(p.coeffs)
Base.lastindex(p::MutableDensePolynomial) = length(p) - 1
function Base.getindex(p::MutableDensePolynomial{B,T,X}, i::Int) where {B,T,X}
    (i < firstindex(p) || i > lastindex(p)) && return zero(T)
    p.coeffs[i + offset(p)]
end
# ??? should this call chop! if `iszero(value)`?
function Base.setindex!(p::P, value, i::Int) where {B,T,X,P<:MutableDensePolynomial{B,T,X}}
    a,b = firstindex(p), lastindex(p)
    iszero(p) && return P([value], i)
    z = zero(first(p.coeffs))
    if i > b
        append!(p.coeffs, _zeros(p, z, i - b))
        p.coeffs[end] = value
    elseif i < a
        prepend!(p.coeffs, _zeros(p, z, a-i))
        p.coeffs[1] = value
    else
        p.coeffs[i + offset(p)] = value
    end
    p
end

offset(p::MutableDensePolynomial) = 1
Base.eachindex(p::MutableDensePolynomial) = 0:1:lastindex(p)
Base.iterate(p::MutableDensePolynomial, args...) = Base.iterate(p.coeffs, args...)
Base.pairs(p::MutableDensePolynomial) =
    Base.Generator(=>, 0:lastindex(p), p.coeffs)

# return a container of zeros based on basis type
_zeros(::Type{<:MutableDensePolynomial}, z, N)  = fill(z, N)

Base.similar(p::MutableDensePolynomial, args...) = similar(p.coeffs, args...)

# iszero
Base.iszero(p::MutableDensePolynomial) = iszero(p.coeffs)::Bool

function degree(p::MutableDensePolynomial)
    length(p.coeffs) - 1
end

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

function chop!(p::MutableDensePolynomial{B,T,X};
               atol=nothing, rtol=Base.rtoldefault(float(real(T)))) where {B,T,X}
    iᵣ = chop_right_index(p.coeffs; atol=atol, rtol=rtol) # nothing -- nothing to chop
    iᵣ === nothing && return zero(p)
    iᵣ == length(p.coeffs) && return p

    N = length(p.coeffs)

    for i ∈ (iᵣ+1):N
        pop!(p.coeffs)
    end
    p
end

function normΔ(q1::MutableDensePolynomial{B}, q2::MutableDensePolynomial{B}) where {B}
    iszero(q1) && return norm(q2,2)
    iszero(q2) && return norm(q1,2)
    r = abs(zero(q1[end] + q2[end]))
    tot = zero(r)
    for i ∈ o:maximum(lastindex, (q1,q2))
       @inbounds tot += abs2(q1[i] - q2[i])
    end
    return sqrt(tot)
end

minimumexponent(::Type{<:MutableDensePolynomial}) = 0



# vector ops +, -, c*x
## unary - (map is as fast)
## binary +
Base.:+(p::MutableDensePolynomial{B,T,X}, q::MutableDensePolynomial{B,S,X}) where{B,X,T,S} =
    _vector_combine(+, p, q)
Base.:-(p::MutableDensePolynomial{B,T,X}, q::MutableDensePolynomial{B,S,X}) where{B,X,T,S} =
    _vector_combine(-, p, q)
function _vector_combine(op, p::MutableDensePolynomial{B,T,X}, q::MutableDensePolynomial{B,S,X}) where{B,X,T,S}
    R = promote_type(T,S)
    P = MutableDensePolynomial{B,R,X}

    iszero(p) && return convert(P, op(q))
    iszero(q) && return convert(P, p)

    b₁, b₂ = lastindex(p), lastindex(q)
    a, b = 0, max(b₁, b₂)

    N = b - a + 1
    z = zero(first(p.coeffs) + first(q.coeffs))
    x = _zeros(p, z, N)

    # zip faster than `pairs`
    @inbounds for (i, cᵢ) ∈ enumerate(p.coeffs)
        x[i] = cᵢ
    end
    @inbounds for (i, cᵢ) ∈ enumerate(q.coeffs)
        x[i] = op(x[i], cᵢ)
    end

    b₁ == b₂ && (x = trim_trailing_zeros(x))
    P(Val(false), x, a)

end
