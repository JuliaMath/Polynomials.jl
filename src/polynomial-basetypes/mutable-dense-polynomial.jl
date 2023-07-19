# * has order
# * leading 0s are trimmed
# * pass Val(true) to bypass trimmings
struct MutableDensePolynomial{B,T,X} <: AbstractUnivariatePolynomial{B,T, X}
    coeffs::Vector{T}
    order::Int # lowest degree, typically 0
    function MutableDensePolynomial{B,T,X}(cs, order::Int=0) where {B,T,X}
        if Base.has_offset_axes(cs)
            @warn "Using the axis offset of the coefficient vector"
            cs, order = cs.parent, first(cs.offsets)
        end

        i = findlast(!iszero, cs)
        if i == nothing
            xs = T[]
        else
            j = findfirst(!iszero, cs)
            xs = T[cs[i] for i ∈ j:i]
            order = order + j - 1
        end
        new{B,T,Symbol(X)}(xs, order)
    end
    function MutableDensePolynomial{B,T,X}(checked::Val{false}, cs::Vector{T}, order::Int=0) where {B,T,X}
        if Base.has_offset_axes(cs)
            @warn "Using the axis offset of the coefficient vector"
            cs, order = cs.parent, first(cs.offsets)
        end
        new{B,T,Symbol(X)}(cs, order)
    end
end


MutableDensePolynomial{B,T,X}(checked::Val{true}, cs::Vector{T}, order::Int=0) where {B,T,X} = MutableDensePolynomial{B,T,X}(cs, order)

function MutableDensePolynomial{B,T}(xs::AbstractVector{S}, order::Int=0, var::SymbolLike=Var(:x)) where {T, S, B}
    if Base.has_offset_axes(xs)
        @warn "Using the axis offset of the coefficient vector"
        xs, order = xs.parent, first(xs.offsets)
    end

    cs = convert(Vector{T}, xs)
    cs = trim_trailing_zeros(cs)
    MutableDensePolynomial{B,T,Symbol(var)}(Val(true),cs, order)
end

function MutableDensePolynomial{B}(xs::AbstractVector{T}, order::Int=0, var::SymbolLike=Var(:x)) where {B, T}
    MutableDensePolynomial{B,T}(xs, order, var)
end

# function MutableDensePolynomial{B,X}(xs::AbstractVector{T}, order::Int) where {T, B,X}
#     MutableDensePolynomial{B,T,X}(xs,order)
# end

# function MutableDensePolynomial{B,X}(xs, order::Int) where {B,X}
#     cs = collect(xs)
#     T = eltype(cs)
#     MutableDensePolynomial{B,T,X}(cs, order)
# end

function MutableDensePolynomial{B,T}(xs, var::SymbolLike) where {B,T}
    cs = convert(Vector{T}, xs)
    MutableDensePolynomial{B,T,Symbol(var)}(cs, 0)
end

# allow specification of codefficients, order, symbol or coefficients, symbol
function MutableDensePolynomial{B}(xs, order::Int=0, var::SymbolLike=Var(:x)) where {B}
    cs = collect(promote(xs...))
    T = eltype(cs)
    MutableDensePolynomial{B, T, Symbol(var)}(cs, order)
end

function MutableDensePolynomial{B}(xs, var::SymbolLike) where {B}
    cs = collect(xs)
    T = eltype(cs)
    MutableDensePolynomial{B, T, Symbol(var)}(cs)
end

function _polynomial(p::P, as::Vector{S})  where {B,T, X, P <: MutableDensePolynomial{B,T,X}, S}
    R = eltype(as)
    Q = MutableDensePolynomial{B, R, X}
    as = trim_trailing_zeros(as)
    Q(Val(false), as, p.order)
end

@poly_register MutableDensePolynomial
constructorof(::Type{<:MutableDensePolynomial{B}}) where {B} = MutableDensePolynomial{B}



## Generics for polynomials
function Base.convert(::Type{MutableDensePolynomial{B,T,X}}, q::MutableDensePolynomial{B,T′,X′}) where {B,T,T′,X,X′}
    MutableDensePolynomial{B,T,X}(Val(false), convert(Vector{T},q.coeffs), q.order)
end

# function Base.convert(::Type{MutableDensePolynomial{B,T}}, q::MutableDensePolynomial{B,T′,X′}) where {B,T,T′,X′}
#     convert(Val(false), MutableDensePolynomial{B,T,X′}, q)
# end


Base.copy(p::MutableDensePolynomial{B,T,X}) where {B,T,X} =
    MutableDensePolynomial{B,T,Var(X)}(copy(p.coeffs), firstindex(p))

# This is B <: StandardBasis?
Base.firstindex(p::MutableDensePolynomial) = p.order
Base.length(p::MutableDensePolynomial) = length(p.coeffs)
offset(p::MutableDensePolynomial) = 1 - firstindex(p)
function Base.getindex(p::MutableDensePolynomial{B,T,X}, i::Int) where {B,T,X}
    (i < firstindex(p) || i > lastindex(p)) && return zero(T)
    p.coeffs[i + offset(p)]
end
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
    P(p.coeffs, o)
end

Base.lastindex(p::MutableDensePolynomial) = firstindex(p) + length(p) - 1
Base.eachindex(p::MutableDensePolynomial) = firstindex(p):1:lastindex(p)
Base.iterate(p::MutableDensePolynomial, args...) = Base.iterate(p.coeffs, args...)
Base.pairs(p::MutableDensePolynomial) =
    Base.Generator(=>, firstindex(p):lastindex(p), p.coeffs)
#    (i - offset(p) => cᵢ for (i, cᵢ) ∈ enumerate(p))
# return a container of zeros based on basis type
_zeros(::Type{<:MutableDensePolynomial}, z, N)  = fill(z, N)

Base.similar(p::MutableDensePolynomial, args...) = similar(p.coeffs, args...)

# iszero, isconstant
Base.iszero(p::MutableDensePolynomial) = iszero(p.coeffs)::Bool

function degree(p::MutableDensePolynomial)
    i = findlast(!iszero, p.coeffs)
    isnothing(i) && return -1
    firstindex(p) + i - 1
end

# zero, one, variable, basis
Base.zero(::Type{MutableDensePolynomial{B,T,X}}) where {B,T,X} =
    MutableDensePolynomial{B,T,X}(T[])

coeffs(p::MutableDensePolynomial) = p.coeffs


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
    iₗ = chop_left_index(p.coeffs; atol=atol, rtol=rtol)
    iₗ === nothing && return zero(p)
    iₗ == 1 && iᵣ == length(p.coeffs) && return p

    N = length(p.coeffs)

    o = order(p)
    for i ∈ 1:(iₗ-1)
        popfirst!(p.coeffs)
        o += 1
    end

    for i ∈ (iᵣ+1):N
        pop!(p.coeffs)
    end

    MutableDensePolynomial{B,T,X}(p.coeffs, o)
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


# vector ops +, -, c*x
## unary
Base.:-(p::MutableDensePolynomial{B,T,X}) where {B,T,X} =
    MutableDensePolynomial{B,T,X}(Val(false), -p.coeffs, p.order)

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

# slower
# function Base.:+(p::MutableDensePolynomial{B,T,X}, q::MutableDensePolynomial{B,S,X}) where {B,T,S,X}
#     a = min(firstindex(p), firstindex(q))
#     b = max(lastindex(p), lastindex(q))
#     x = collect(p[i] - q[i] for i ∈ a:b)
#     MutableDensePolynomial{B,eltype(x),X}(x,a)
# end

# scalar
function scalar_mult(p::MutableDensePolynomial{B,T,X}, c::S) where {B,T,X,S}
    cs = p.coeffs .* (c,) # works with T[]
    return _polynomial(p, cs)
end
function scalar_mult(c::S, p::MutableDensePolynomial{B,T,X}) where {B,T,X,S}
    cs = (c,) .* p.coeffs
    return _polynomial(p, cs)
end

#XXX need to add LinearAlgebra
# function LinearAlgebra.lmul!(c::Number, p::MutableDensePolynomial{B,T,X}) where {B,T,X}
#     MutableDensePolynomial{B,T,X}(Val(false), lmul!(c, p.coeffs), firstindex(p))
# end

# function LinearAlgebra.rmul!(c::Number, p::MutableDensePolynomial{B,T,X}) where {B,T,X}
#     MutableDensePolynomial{B,T,X}(rmul!(c, p.coeffs), firstindex(p))
# end
