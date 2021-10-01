using LinearAlgebra

export fromroots,
       truncate!,
       chop!,
       coeffs,
       degree,
       domain,
       mapdomain,
       order,
       hasnan,
       roots,
       companion,
       vander,
       fit,
       integrate,
       derivative,
       variable,
       isintegral,
       ismonic

"""
    fromroots(::AbstractVector{<:Number}; var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractVector{<:Number}; var=:x)

Construct a polynomial of the given type given the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest common
julia> using Polynomials

julia> r = [3, 2]; # (x - 3)(x - 2)

julia> fromroots(r)
Polynomial(6 - 5*x + x^2)
```
"""
function fromroots(P::Type{<:AbstractPolynomial}, roots::AbstractVector; var::SymbolLike = :x)
    x = variable(P, var)
    p =  prod(x - r for r in roots)
    return truncate!(p)
end
fromroots(r::AbstractVector{<:Number}; var::SymbolLike = :x) =
    fromroots(Polynomial, r, var = var)

"""
    fromroots(::AbstractMatrix{<:Number}; var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractMatrix{<:Number}; var=:x)

Construct a polynomial of the given type using the eigenvalues of the given matrix as the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest common
julia> using Polynomials

julia> A = [1 2; 3 4]; # (x - 5.37228)(x + 0.37228)

julia> fromroots(A)
Polynomial(-1.9999999999999998 - 5.0*x + 1.0*x^2)
```
"""
fromroots(P::Type{<:AbstractPolynomial},
    A::AbstractMatrix{T};
    var::SymbolLike = :x,) where {T <: Number} = fromroots(P, eigvals(A), var = var)
fromroots(A::AbstractMatrix{T}; var::SymbolLike = :x) where {T <: Number} =
    fromroots(Polynomial, eigvals(A), var = var)

"""
    fit(x, y, deg=length(x) - 1; [weights], var=:x)
    fit(::Type{<:AbstractPolynomial}, x, y, deg=length(x)-1; [weights], var=:x)

Fit the given data as a polynomial type with the given degree. Uses
linear least squares to minimize the norm  `||y - V⋅β||^2`, where `V` is
the Vandermonde matrix and `β` are the coefficients of the polynomial
fit.

This will automatically scale your data to the [`domain`](@ref) of the
polynomial type using [`mapdomain`](@ref). The default polynomial type
is [`Polynomial`](@ref).


## Weights

Weights may be assigned to the points by specifying a vector or matrix of weights.

When specified as a vector, `[w₁,…,wₙ]`, the weights should be
non-negative as the minimization problem is `argmin_β Σᵢ wᵢ |yᵢ - Σⱼ
Vᵢⱼ βⱼ|² = argmin_β || √(W)⋅(y - V(x)β)||²`, where, `W` the digonal
matrix formed from `[w₁,…,wₙ]`, is used for the solution, `V` being
the Vandermonde matrix of `x` corresponding to the specified
degree. This parameterization of the weights is different from that of
`numpy.polyfit`, where the weights would be specified through
`[ω₁,ω₂,…,ωₙ] = [√w₁, √w₂,…,√wₙ]`
with the answer solving
`argminᵦ | (ωᵢ⋅yᵢ- ΣⱼVᵢⱼ(ω⋅x) βⱼ) |^2`.

When specified as a matrix, `W`, the solution is through the normal
equations `(VᵀWV)β = (Vᵀy)`, again `V` being the Vandermonde matrix of
`x` corresponding to the specified degree.

(In statistics, the vector case corresponds to weighted least squares,
where weights are typically given by `wᵢ = 1/σᵢ²`, the `σᵢ²` being the
variance of the measurement; the matrix specification follows that of
the generalized least squares estimator with `W = Σ⁻¹`, the inverse of
the variance-covariance matrix.)

## large degree

For fitting with a large degree, the Vandermonde matrix is exponentially ill-conditioned. The [`ArnoldiFit`](@ref) type introduces an Arnoldi orthogonalization that fixes this problem.

"""
function fit(P::Type{<:AbstractPolynomial},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::Integer = length(x) - 1;
             weights = nothing,
             var = :x,) where {T}
    _fit(P, x, y, deg; weights=weights, var=var)
end

fit(P::Type{<:AbstractPolynomial},
    x,
    y,
    deg::Integer = length(x) - 1;
    weights = nothing,
    var = :x,) = fit′(P, promote(collect(x), collect(y))..., deg; weights = weights, var = var)

#  avoid issue  214
fit′(P::Type{<:AbstractPolynomial}, x, y, args...;kwargs...) = throw(ArgumentError("x and y do not produce abstract   vectors"))
fit′(P::Type{<:AbstractPolynomial},
     x::AbstractVector{T},
     y::AbstractVector{T},
     args...; kwargs...) where {T} = fit(P,x,y,args...;  kwargs...)


fit(x::AbstractVector,
    y::AbstractVector,
    deg::Integer = length(x) - 1;
    weights = nothing,
    var = :x,) = fit(Polynomial, x, y, deg; weights = weights, var = var)

function _fit(P::Type{<:AbstractPolynomial},
             x::AbstractVector{T},
             y::AbstractVector{T},
             deg::Integer = length(x) - 1;
             weights = nothing,
             var = :x,) where {T}
    x = mapdomain(P, x)
    vand = vander(P, x, deg)
    if weights !== nothing
        coeffs = _wlstsq(vand, y, weights)
    else
        coeffs = qr(vand) \ y
    end
    R = float(T)
    return P(R.(coeffs), var)
end


# Weighted linear least squares
_wlstsq(vand, y, W::Number) = _wlstsq(vand, y, fill!(similar(y), W))
function _wlstsq(vand, y, w::AbstractVector)
    W = Diagonal(sqrt.(w))
    qr(W * vand) \ (W * y)
end
_wlstsq(vand, y, W::AbstractMatrix) = qr(vand' * W * vand) \ (vand' * W * y)

"""
    roots(::AbstractPolynomial; kwargs...)

Returns the roots, or zeros, of the given polynomial.

For non-factored, standard basis polynomials the roots are calculated via the eigenvalues of the companion matrix. The `kwargs` are passed to the `LinearAlgeebra.eigvals` call.

!!! note
    The default `roots` implementation is for polynomials in the
    standard basis. The companion matrix approach is reasonably fast
    and accurate for modest-size polynomials. However, other packages
    in the `Julia` ecosystem may be of interest and are mentioned in the documentation.


"""
function roots(q::AbstractPolynomial{T}; kwargs...) where {T}

    p = convert(Polynomial{T},  q)
    roots(p; kwargs...)

end

"""
    companion(::AbstractPolynomial)

Return the companion matrix for the given polynomial.

# References
[Companion Matrix](https://en.wikipedia.org/wiki/Companion_matrix)
"""
companion(::AbstractPolynomial)

"""
    vander(::Type{AbstractPolynomial}, x::AbstractVector, deg::Integer)

Calculate the pseudo-Vandermonde matrix of the given polynomial type with the given degree.

# References
[Vandermonde Matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix)
"""
vander(::Type{<:AbstractPolynomial}, x::AbstractVector, deg::Integer)

"""
    integrate(p::AbstractPolynomial)

Return an antiderivative for `p`
"""
integrate(P::AbstractPolynomial) = throw(ArgumentError("`integrate` not implemented for polynomials of type $P"))

"""
    integrate(::AbstractPolynomial, C)

Returns the indefinite integral of the polynomial with constant `C` when expressed in the standard basis.
"""
function integrate(p::P, C) where {P <: AbstractPolynomial}
    ∫p = integrate(p)
    isnan(C) && return ⟒(P){eltype(∫p+C), indeterminate(∫p)}([C])
    ∫p + (C - constantterm(∫p))
end

"""
    integrate(::AbstractPolynomial, a, b)

Compute the definite integral of the given polynomial from `a` to `b`. Will throw an error if either `a` or `b` are out of the polynomial's domain.
"""
function integrate(p::AbstractPolynomial, a, b)
    P = integrate(p)
    return P(b) - P(a)
end

"""
    derivative(::AbstractPolynomial, order::Int = 1)

Returns a polynomial that is the `order`th derivative of the given polynomial. `order` must be non-negative.
"""
derivative(::AbstractPolynomial, ::Int)

"""
    truncate!(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

In-place version of [`truncate`](@ref)
"""
function truncate!(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {T}
    truncate!(p.coeffs, rtol=rtol, atol=atol)
    return chop!(p, rtol = rtol, atol = atol)
end

## truncate! underlying storage type
function truncate!(ps::Vector{T};
                   rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {T}
    max_coeff = norm(ps, Inf)
    thresh = max_coeff * rtol + atol
    for (i,pᵢ) ∈ pairs(ps)
        if abs(pᵢ) <= thresh
            ps[i] = zero(T)
        end
    end
    nothing
end

function truncate!(ps::Dict{Int,T};
                   rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {T}

    max_coeff = norm(values(ps), Inf)
    thresh = max_coeff * rtol + atol

    for (k,val) in  ps
        if abs(val) <= thresh
            pop!(ps,k)
        end
    end
    nothing
end

truncate!(ps::NTuple; kwargs...) = throw(ArgumentError("`truncate!` not defined."))
Base.truncate(ps::NTuple{0}; kwargs...) = ps
function Base.truncate(ps::NTuple{N,T};
              rtol::Real = Base.rtoldefault(real(T)),
              atol::Real = 0,) where {N,T}
    thresh = norm(ps, Inf) * rtol + atol
    return NTuple{N,T}(abs(pᵢ) <= thresh ? zero(T) : pᵢ for pᵢ ∈ values(ps))
end


"""
    truncate(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

Rounds off coefficients close to zero, as determined by `rtol` and `atol`, and then chops any leading zeros. Returns a new polynomial.
"""
function Base.truncate(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
    atol::Real = 0,) where {T}
    truncate!(deepcopy(p), rtol = rtol, atol = atol)
end

"""
    chop!(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0))

In-place version of [`chop`](@ref)
"""
function chop!(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
               atol::Real = 0,) where {T}
    chop!(p.coeffs, rtol=rtol, atol=atol)
    return p
end

# chop! underlying storage type
function chop!(ps::Vector{T};
               rtol::Real = Base.rtoldefault(real(T)),
               atol::Real = 0,) where {T}

    tol = norm(ps) * rtol + atol
    for i = lastindex(ps):-1:1
        val = ps[i]
        if abs(val) > tol
            resize!(ps, i);
            return nothing
        end
    end
    resize!(ps, 1)
    return nothing
end

function chop!(ps::Dict{Int,T};
    rtol::Real = Base.rtoldefault(real(T)),
    atol::Real = 0,) where {T}

    tol = norm(values(ps)) * rtol + atol

    for k in sort(collect(keys(ps)), by=x->x[1], rev=true)
        if  abs(ps[k]) > tol
            return nothing
        end
        pop!(ps, k)
    end

    return nothing
end

chop!(ps::NTuple; kwargs...) = throw(ArgumentError("chop! not defined"))
Base.chop(ps::NTuple{0}; kwargs...) = ps
function Base.chop(ps::NTuple{N,T};
              rtol::Real = Base.rtoldefault(real(T)),
              atol::Real = 0,) where {N,T}
    thresh = norm(ps, Inf) * rtol + atol
    for i in N:-1:1
        if abs(ps[i]) > thresh
            return ps[1:i]
        end
    end
    return NTuple{0,T}()
end



"""
    chop(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0))

Removes any leading coefficients that are approximately 0 (using `rtol` and `atol`). Returns a polynomial whose degree will guaranteed to be equal to or less than the given polynomial's.
"""
function Base.chop(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
    atol::Real = 0,) where {T}
    chop!(deepcopy(p), rtol = rtol, atol = atol)
end




"""
    check_same_variable(p::AbstractPolynomial, q::AbstractPolynomial)

Check if either `p` or `q` is constant or if `p` and `q` share the same variable
"""
check_same_variable(p::AbstractPolynomial, q::AbstractPolynomial) =
    (isconstant(p) || isconstant(q)) || indeterminate(p) ==  indeterminate(q)

function assert_same_variable(p::AbstractPolynomial, q::AbstractPolynomial)
    check_same_variable(p,q) || throw(ArgumentError("Polynomials have different indeterminates"))
end

function assert_same_variable(X::Symbol, Y::Symbol)
    X == Y ||  throw(ArgumentError("Polynomials have different indeterminates"))
end

#=
Linear Algebra =#
"""
    norm(::AbstractPolynomial, p=2)

Calculates the p-norm of the polynomial's coefficients
"""
function LinearAlgebra.norm(q::AbstractPolynomial, p::Real = 2)
    vs = values(q)
    isempty(vs) && return zero(eltype(q))
    return norm(vs, p) # if vs=() must be handled in special type
end

"""
    conj(::AbstractPolynomial)

Returns the complex conjugate of the polynomial
"""
LinearAlgebra.conj(p::P) where {P <: AbstractPolynomial} = map(conj, p)
LinearAlgebra.adjoint(p::P) where {P <: AbstractPolynomial} = map(adjoint, p)
LinearAlgebra.adjoint(A::VecOrMat{<:AbstractPolynomial}) = adjoint.(permutedims(A)) # default has indeterminate indeterminate
LinearAlgebra.transpose(p::AbstractPolynomial) = p
LinearAlgebra.transpose!(p::AbstractPolynomial) = p

#=
Conversions =#
Base.convert(::Type{P}, p::P) where {P <: AbstractPolynomial} = p
Base.convert(P::Type{<:AbstractPolynomial}, x) = P(x)
function Base.convert(::Type{S}, p::P) where {S <: Number,T, P<:AbstractPolynomial{T}}
    isconstant(p) && return convert(S, constantterm(p))
    throw(ArgumentError("Can't convert a nonconstant polynomial to type $S"))
end
function Base.convert(::Type{T}, p::P) where {T, P<:AbstractPolynomial{T}}
    isconstant(p) && return constantterm(p)
    throw(ArgumentError("Can't convert a nonconstant polynomial to type $S"))
end

# Methods to ensure that matrices of polynomials behave as desired
Base.promote_rule(::Type{<:AbstractPolynomial{T}},
                  ::Type{<:AbstractPolynomial{S}},
                  ) where {T,S} = Polynomial{promote_type(T, S)}
Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,   Q<:AbstractPolynomial{S,X}} =
                                                   Polynomial{promote_type(T, S),X}
Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,Y, Q<:AbstractPolynomial{S,Y}} =
                                                  assert_same_variable(X,Y)


#=
Inspection =#
"""
    length(::AbstractPolynomial)

The length of the polynomial.
"""
Base.length(p::AbstractPolynomial) = length(coeffs(p))
"""
    size(::AbstractPolynomial, [i])

Returns the size of the polynomials coefficients, along axis `i` if provided.
"""
Base.size(p::AbstractPolynomial) = size(coeffs(p))
Base.size(p::AbstractPolynomial, i::Integer) = size(coeffs(p), i)
Base.eltype(p::AbstractPolynomial{T}) where {T} = T
# in  analogy  with  polynomial as a Vector{T} with different operations defined.
Base.eltype(::Type{<:AbstractPolynomial}) = Float64
Base.eltype(::Type{<:AbstractPolynomial{T}}) where {T} = T
_eltype(::Type{<:AbstractPolynomial}) = nothing
_eltype(::Type{<:AbstractPolynomial{T}}) where {T} = T
function _eltype(P::Type{<:AbstractPolynomial}, p::AbstractPolynomial)
    T′ = _eltype(P)
    T = T′ == nothing ? eltype(p) : T′
    T
end
Base.iszero(p::AbstractPolynomial) = all(iszero, p)


# See discussions in https://github.com/JuliaMath/Polynomials.jl/issues/258
"""
    all(pred, poly::AbstractPolynomial)

Test whether all coefficients of an `AbstractPolynomial` satisfy predicate `pred`.

You can implement `isreal`, etc., to a `Polynomial` by using `all`.
"""
Base.all(pred, p::AbstractPolynomial) = all(pred, values(p))
"""
    any(pred, poly::AbstractPolynomial)

Test whether any coefficient of an `AbstractPolynomial` satisfies predicate `pred`.
"""
Base.any(pred, p::AbstractPolynomial) = any(pred, values(p))




"""
    map(fn, p::AbstractPolynomial, args...)

Transform coefficients of `p` by applying a function (or other callables) `fn` to each of them.

You can implement `real`, etc., to a `Polynomial` by using `map`.
"""
Base.map(fn, p::P, args...)  where {P<:AbstractPolynomial} = _convert(p, map(fn, coeffs(p), args...))

"""
    isreal(p::AbstractPolynomial)

Determine whether a polynomial is a real polynomial, i.e., having only real numbers as coefficients.

See also: [`real`](@ref)
"""
Base.isreal(p::AbstractPolynomial) = all(isreal, p)
"""
    real(p::AbstractPolynomial)

Construct a real polynomial from the real parts of the coefficients of `p`.

See also: [`isreal`](@ref)

!!! note
    This could cause losing terms in `p`. This method is usually called on polynomials like `p = Polynomial([1, 2 + 0im, 3.0, 4.0 + 0.0im])` where you want to chop the imaginary parts of the coefficients of `p`.
"""
Base.real(p::AbstractPolynomial) = map(real, p)

"""
    isintegral(p::AbstractPolynomial)

Determine whether a polynomial is an integer polynomial, i.e., having only integers as coefficients.
"""
isintegral(p::AbstractPolynomial) = all(isinteger, p)

"""
    ismonic(p::AbstractPolynomial)

Determine whether a polynomial is a monic polynomial, i.e., its leading coefficient is one.
"""
ismonic(p::AbstractPolynomial) = isone(convert(Polynomial, p)[end])

"`hasnan(p::AbstractPolynomial)` are any coefficients `NaN`"
hasnan(p::AbstractPolynomial) = any(hasnan, p)
hasnan(p::AbstractArray) = any(hasnan.(p))
hasnan(x) = isnan(x)

"""
    isconstant(::AbstractPolynomial)

Is the polynomial  `p` a constant.
"""
isconstant(p::AbstractPolynomial) = degree(p) <= 0

"""
    coeffs(::AbstractPolynomial)

Return the coefficient vector. For a standard basis polynomial these are `[a_0, a_1, ..., a_n]`.
"""
coeffs(p::AbstractPolynomial) = p.coeffs

# hook in for offset coefficients of Laurent Polynomials
_coeffs(p::AbstractPolynomial) = coeffs(p)


# specialize this to p[0] when basis vector is 1
"""
    constantterm(p::AbstractPolynomial)

return `p(0)`, the constant term in the standard basis
"""
constantterm(p::AbstractPolynomial{T}) where {T} = p(zero(T))

"""
    degree(::AbstractPolynomial)

Return the degree of the polynomial, i.e. the highest exponent in the polynomial that
has a nonzero coefficient. The degree of the zero polynomial is defined to be -1. The default method assumes the basis polynomial, `βₖ` has degree `k`.
"""
degree(p::AbstractPolynomial) = iszero(p) ? -1 : lastindex(p)


"""
    domain(::Type{<:AbstractPolynomial})

Returns the domain of the polynomial.
"""
domain(::Type{<:AbstractPolynomial})
domain(::P) where {P <: AbstractPolynomial} = domain(P)

"""
    mapdomain(::Type{<:AbstractPolynomial}, x::AbstractArray)
    mapdomain(::AbstractPolynomial, x::AbstractArray)

Given values of x that are assumed to be unbounded (-∞, ∞), return values rescaled to the domain of the given polynomial.

# Examples
```jldoctest  common
julia> using Polynomials

julia> x = -10:10
-10:10

julia> extrema(mapdomain(ChebyshevT, x))
(-1.0, 1.0)

```
"""
function mapdomain(P::Type{<:AbstractPolynomial}, x::AbstractArray)
    d = domain(P)
    x = collect(x)
    x_zerod = x .- minimum(x)
    x_scaled = x_zerod .* (last(d) - first(d)) ./ maximum(x_zerod)
    x_scaled .+= first(d)
    return x_scaled
end
mapdomain(::P, x::AbstractArray) where {P <: AbstractPolynomial} = mapdomain(P, x)

#=
indexing =#
Base.firstindex(p::AbstractPolynomial) = 0
Base.lastindex(p::AbstractPolynomial) = length(p) - 1
Base.eachindex(p::AbstractPolynomial) = 0:length(p) - 1
Base.broadcastable(p::AbstractPolynomial) = Ref(p)

# getindex
function Base.getindex(p::AbstractPolynomial{T}, idx::Int) where {T}
    m,M = firstindex(p), lastindex(p)
    idx < m && throw(BoundsError(p, idx))
    idx > M && return zero(T)
    p.coeffs[idx - m + 1]
end
Base.getindex(p::AbstractPolynomial, idx::Number) = getindex(p, convert(Int, idx))
Base.getindex(p::AbstractPolynomial, indices) = [p[i] for i in indices]
Base.getindex(p::AbstractPolynomial, ::Colon) = coeffs(p)

# setindex
function Base.setindex!(p::AbstractPolynomial, value, idx::Int)
    n = length(coeffs(p))
    if n ≤ idx
        resize!(p.coeffs, idx + 1)
        p.coeffs[n + 1:idx] .= 0
    end
    p.coeffs[idx + 1] = value
    return p
end

Base.setindex!(p::AbstractPolynomial, value, idx::Number) =
    setindex!(p, value, convert(Int, idx))
#Base.setindex!(p::AbstractPolynomial, value::Number, indices) =
#    [setindex!(p, value, i) for i in indices]
Base.setindex!(p::AbstractPolynomial, values, indices) =
    [setindex!(p, v, i) for (v, i) in tuple.(values, indices)]
#    [setindex!(p, v, i) for (v, i) in zip(values, indices)]
#Base.setindex!(p::AbstractPolynomial, value, ::Colon) =
#    setindex!(p, value, eachindex(p))
Base.setindex!(p::AbstractPolynomial, values, ::Colon) =
#        [setindex!(p, v, i) for (v, i) in zip(values, eachindex(p))]
    [setindex!(p, v, i) for (v, i) in tuple.(values, eachindex(p))]

#=
Iteration =#
## XXX breaking change in v2.0.0
#
# we want to keep iteration close to iteration over the coefficients as a vector:
# `iterate` iterates over coefficients, includes 0s
# `collect(T, p)` yields coefficients of `p` converted to type `T`
# `keys(p)` an iterator spanning `firstindex` to `lastindex` which *may* skip 0 terms (SparsePolynomial)
#    and *may* not be in order (SparsePolynomial)
# `values(p)` `pᵢ` for each `i` in `keys(p)`
# `pairs(p)`: `i => pᵢ` possibly skipping over values of `i` with `pᵢ == 0` (SparsePolynomial)
#    and possibly non ordered (SparsePolynomial)
# `monomials(p)`: iterates over pᵢ ⋅ basis(p, i) i  ∈ keys(p)
function Base.iterate(p::AbstractPolynomial, state=nothing)
    i = firstindex(p)
    if state == nothing
        return (p[i], i)
    else
        j = lastindex(p)
        if i <= state < j
            return (p[state+1], state+1)
        end
        return nothing
    end
end

# pairs map i -> aᵢ *possibly* skipping over ai == 0
# cf. abstractdict.jl
struct PolynomialKeys{P}
    p::P
end
struct PolynomialValues{P}
    p::P
end
Base.keys(p::AbstractPolynomial) =  PolynomialKeys(p)
Base.values(p::AbstractPolynomial) =  PolynomialValues(p)
Base.length(p::PolynomialValues) = length(p.p.coeffs)
Base.length(p::PolynomialKeys) = length(p.p.coeffs)
Base.size(p::Union{PolynomialValues, PolynomialKeys}) = (length(p),)
function Base.iterate(v::PolynomialKeys, state=nothing)
    i = firstindex(v.p)
    state==nothing && return (i, i)
    j = lastindex(v.p)
    i <= state < j && return (state+1, state+1)
    return nothing
end

function Base.iterate(v::PolynomialValues, state=nothing)
    i = firstindex(v.p)
    state==nothing && return (v.p[i], i)
    j = lastindex(v.p)
    i <= state < j && return (v.p[state+1], state+1)
    return nothing
end


# iterate over monomials of the polynomial
struct Monomials{P}
    p::P
end
"""
    monomials(p::AbstractPolynomial)

Returns an iterator over the terms, `pᵢ⋅basis(p,i)`, of the polynomial for each `i` in `keys(p)`.
"""
monomials(p) = Monomials(p)
function Base.iterate(v::Monomials, state...)
    y = iterate(pairs(v.p), state...)
    y == nothing && return nothing
    kv, s = y
    return (kv[2]*basis(v.p, kv[1]), s)
end
Base.length(v::Monomials) = length(keys(v.p))


#=
identity =#
Base.copy(p::P) where {P <: AbstractPolynomial} = _convert(p, copy(coeffs(p)))
Base.hash(p::AbstractPolynomial, h::UInt) = hash(indeterminate(p), hash(coeffs(p), h))

# get symbol of polynomial. (e.g. `:x` from 1x^2 + 2x^3...
#_indeterminate(::Type{P}) where {T, X, P <: AbstractPolynomial{T, X}} = X
_indeterminate(::Type{P}) where {P <: AbstractPolynomial} = nothing
_indeterminate(::Type{P}) where {T, X, P <: AbstractPolynomial{T,X}} = X
function indeterminate(::Type{P}) where {P <: AbstractPolynomial}
    X = _indeterminate(P)
    X == nothing ? :x : X
end
indeterminate(p::P) where {P <: AbstractPolynomial} = _indeterminate(P)
function indeterminate(PP::Type{P}, p::AbstractPolynomial{T,Y}) where {P <: AbstractPolynomial, T,Y}
    X = _indeterminate(PP)
    X == nothing && return Y
    assert_same_variable(X,Y)
    return X
    #X = _indeterminate(PP) == nothing ? indeterminate(p) :  _indeterminate(PP)
end
function indeterminate(PP::Type{P}, x::Symbol) where {P <: AbstractPolynomial}
    X = _indeterminate(PP) == nothing ? x :  _indeterminate(PP)
end

#=
zero, one, variable, basis =#



"""
    zero(::Type{<:AbstractPolynomial})
    zero(::AbstractPolynomial)

Returns a representation of 0 as the given polynomial.
"""
function Base.zero(::Type{P}) where {P<:AbstractPolynomial}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}(zeros(T,1))
end
Base.zero(::Type{P}, var::SymbolLike) where {P <: AbstractPolynomial} = zero(⟒(P){eltype(P),Symbol(var)}) #default 0⋅b₀
Base.zero(p::P, var=indeterminate(p)) where {P <: AbstractPolynomial} = zero(P, var)
"""
    one(::Type{<:AbstractPolynomial})
    one(::AbstractPolynomial)

Returns a representation of 1 as the given polynomial.
"""
Base.one(::Type{P}) where {P<:AbstractPolynomial} = throw(ArgumentError("No default method defined")) # no default method
Base.one(::Type{P}, var::SymbolLike) where {P <: AbstractPolynomial} = one(⟒(P){eltype(P), Symbol(var == nothing ? :x : var)})
Base.one(p::P, var=indeterminate(p)) where {P <: AbstractPolynomial} = one(P, var)

Base.oneunit(::Type{P}, args...) where {P <: AbstractPolynomial} = one(P, args...)
Base.oneunit(p::P, args...) where {P <: AbstractPolynomial} = one(p, args...)


"""
    variable(var=:x)
    variable(::Type{<:AbstractPolynomial}, var=:x)
    variable(p::AbstractPolynomial, var=indeterminate(p))

Return the monomial `x` in the indicated polynomial basis.  If no type is give, will default to [`Polynomial`](@ref). Equivalent  to  `P(var)`.

# Examples
```jldoctest  common
julia> using Polynomials

julia> x = variable()
Polynomial(x)

julia> p = 100 + 24x - 3x^2
Polynomial(100 + 24*x - 3*x^2)

julia> roots((x - 3) * (x + 2))
2-element Vector{Float64}:
 -2.0
  3.0

```
"""
variable(::Type{P}) where {P <: AbstractPolynomial} = throw(ArgumentError("No default method defined")) # no default
variable(::Type{P}, var::SymbolLike) where {P <: AbstractPolynomial} = variable(⟒(P){eltype(P),Symbol(var)})
variable(p::AbstractPolynomial, var = indeterminate(p)) = variable(typeof(p), var)
variable(var::SymbolLike = :x) = variable(Polynomial{Int}, var)

# basis
# var is a positional argument, not a keyword; can't deprecate so we do `_var; var=_var`
# return the kth basis polynomial for the given polynomial type, e.g. x^k for Polynomial{T}
function basis(::Type{P}, k::Int) where {P<:AbstractPolynomial}
    T,X = eltype(P), indeterminate(P)
    zs = zeros(T, k+1)
    zs[end] = one(T)
    ⟒(P){eltype(P), indeterminate(P)}(zs)
end
function basis(::Type{P}, k::Int, _var::SymbolLike; var=_var) where {P <: AbstractPolynomial}
    T,X = eltype(P), Symbol(var)
    basis(⟒(P){T,X}, k)
end
basis(p::P, k::Int, _var=indeterminate(p); var=_var) where {P<:AbstractPolynomial} = basis(P, k, var)

#=
arithmetic =#
Base.:-(p::P) where {P <: AbstractPolynomial} = _convert(p, -coeffs(p))

Base.:*(p::AbstractPolynomial, c::Number) = scalar_mult(p, c)
Base.:*(c::Number, p::AbstractPolynomial) = scalar_mult(c, p)
Base.:*(c::T, p::P) where {T, P <: AbstractPolynomial{T}} = scalar_mult(c, p)
Base.:*(p::P, c::T) where {T, P <: AbstractPolynomial{T}} = scalar_mult(p, c)

# implicitly identify c::Number with a constant polynomials
Base.:+(c::Number, p::AbstractPolynomial) = +(p, c)
Base.:-(p::AbstractPolynomial, c::Number) = +(p, -c)
Base.:-(c::Number, p::AbstractPolynomial) = +(-p, c)

# scalar operations
# no generic p+c, as polynomial addition falls back to scalar ops
#function Base.:+(p::P, n::Number) where {P <: AbstractPolynomial}
#    p1, p2 = promote(p, n)
#    return p1 + p2
#end


Base.:-(p1::AbstractPolynomial, p2::AbstractPolynomial) = +(p1, -p2)


## addition
## Fall back addition is possible as vector addition with padding by 0s
## Subtypes will likely want to implement both:
## +(p::P,c::Number) and +(p::P, q::Q) where {T,S,X,P<:SubtypePolynomial{T,X},Q<:SubtypePolynomial{S,X}}
## though the default for poly+poly isn't terrible

# polynomial + scalar; implicit identification of c with c*one(P)
Base.:+(p::P, c::T) where {T,X, P<:AbstractPolynomial{T,X}} = p + c * one(P)

function Base.:+(p::P, c::S) where {T,X, P<:AbstractPolynomial{T,X}, S}
    R = promote_type(T,S)
    q = convert(⟒(P){R,X}, p)
    q + R(c)
end

# polynomial + polynomial when different types
function Base.:+(p::P, q::Q) where {T,X,P <: AbstractPolynomial{T,X}, S,Y,Q <: AbstractPolynomial{S,Y}}
    isconstant(p) && return constantterm(p) + q
    isconstant(q) && return p + constantterm(q)
    assert_same_variable(X,Y)
    sum(promote(p,q))

end

# Works when p,q of same type.
# For Immutable, must remove N,M bit;
# for others, can widen to Type{T,X}, Type{S,X} to avoid a promotion
function Base.:+(p::P, q::P) where {T,X,P<:AbstractPolynomial{T,X}}
    cs = degree(p) >= degree(q)  ? ⊕(P, p.coeffs, q.coeffs) : ⊕(P, q.coeffs, p.coeffs)
    return P(cs)
end

# addition of polynomials is just vector space addition, so can be done regardless
# of basis, as long as the same. These ⊕ methods try to find a performant means to add
# to sets of coefficients based on the storage type. These assume n1 >= n2
function ⊕(P::Type{<:AbstractPolynomial}, p1::Vector{T}, p2::Vector{S}) where {T,S}

    n1, n2 = length(p1), length(p2)
    R = promote_type(T,S)

    cs = collect(R,p1)
    for i in 1:n2
        cs[i] += p2[i]
    end

    return cs
end

# Padded vector sum of two tuples assuming N ≥ M
@generated function ⊕(P::Type{<:AbstractPolynomial}, p1::NTuple{N,T}, p2::NTuple{M,S}) where {T,N,S,M}

    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(p1[$i] + p2[$i])
    end
    for i in M+1:N
        exprs[i] =:(p1[$i])
    end

    return quote
        Base.@_inline_meta
        tuple($(exprs...))
    end

end

# addition when a dictionary is used for storage
function ⊕(P::Type{<:AbstractPolynomial}, p1::Dict{Int,T}, p2::Dict{Int,S}) where {T,S}

    R = promote_type(T,S)
    p = Dict{Int, R}()


    # this allocates in the union
#    for i in union(eachindex(p1), eachindex(p2))
#        p[i] = p1[i] + p2[i]
#    end

    for (i,pi) ∈ pairs(p1)
        @inbounds p[i] = pi + get(p2, i, zero(R))
    end
    for (i,pi) ∈ pairs(p2)
        if iszero(get(p,i,zero(R)))
            @inbounds p[i] = get(p1, i, zero(R)) + pi
        end
    end

    return  p

end

## -- multiplication


# this fall back not necessarily efficient (e.g., sparse)
function scalar_mult(p::P, c::S) where {S, T, X, P<:AbstractPolynomial{T,X}}
    R = Base.promote_op(*, T, S) # typeof(one(T)*one(S))?
    𝐏 = ⟒(P){R,X}
    𝐏([pᵢ * c for pᵢ ∈ coeffs(p)])
end

function scalar_mult(c::S, p::P) where {S, T, X, P<:AbstractPolynomial{T, X}}
    R = Base.promote_op(*, T, S)
    𝐏 = ⟒(P){R,X}
    𝐏([c * pᵢ for pᵢ ∈ coeffs(p)])
end


function Base.:/(p::P, c::S) where {P <: AbstractPolynomial,S}
    _convert(p, coeffs(p) ./ c)
end

## polynomial p*q
## Polynomial multiplication formula depend on the particular basis used. The subtype must implement
function Base.:*(p1::P, p2::Q) where {T,X,P <: AbstractPolynomial{T,X},S,Y,Q <: AbstractPolynomial{S,Y}}
    isconstant(p1) && return constantterm(p1) * p2
    isconstant(p2) && return p1 * constantterm(p2)
    assert_same_variable(X, Y)
    p1, p2 = promote(p1, p2)
    return p1 * p2
end

Base.:^(p::AbstractPolynomial, n::Integer) = Base.power_by_squaring(p, n)

function Base.divrem(num::P, den::O) where {P <: AbstractPolynomial,O <: AbstractPolynomial}
    n, d = promote(num, den)
    return divrem(n, d)
end

"""
    gcd(a::AbstractPolynomial, b::AbstractPolynomial; atol::Real=0, rtol::Real=Base.rtoldefault)

Find the greatest common denominator of two polynomials recursively using
[Euclid's algorithm](http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm).

# Examples

```jldoctest common
julia> using Polynomials

julia> gcd(fromroots([1, 1, 2]), fromroots([1, 2, 3]))
Polynomial(4.0 - 6.0*x + 2.0*x^2)

```
"""
function Base.gcd(p1::AbstractPolynomial{T}, p2::AbstractPolynomial{S}; kwargs...) where {T,S}
    gcd(promote(p1, p2)...; kwargs...)
end

function Base.gcd(p1::AbstractPolynomial{T}, p2::AbstractPolynomial{T};
                  atol::Real=zero(real(T)),
                  rtol::Real=Base.rtoldefault(real(T))
                  ) where {T}


    r₀, r₁ = p1, p2
    iter = 1
    itermax = length(r₁)

    while !iszero(r₁) && iter ≤ itermax
        _, rtemp = divrem(r₀, r₁)
        r₀ = r₁
        r₁ = truncate(rtemp; atol=atol, rtol=rtol)
        iter += 1
    end
    return r₀
end

"""
    uvw(p,q; kwargs...)

return `u` the gcd of `p` and `q`, and `v` and `w`, where `u*v = p` and `u*w = q`.
"""
uvw(p::AbstractPolynomial, q::AbstractPolynomial; kwargs...) = uvw(promote(p,q)...; kwargs...)
uvw(p1::P, p2::P; kwargs...) where {P <:AbstractPolynomial} = throw(ArgumentError("uvw not defind"))

"""
    div(::AbstractPolynomial, ::AbstractPolynomial)
"""
Base.div(n::AbstractPolynomial, d::AbstractPolynomial) = divrem(n, d)[1]

"""
    rem(::AbstractPolynomial, ::AbstractPolynomial)
"""
Base.rem(n::AbstractPolynomial, d::AbstractPolynomial) = divrem(n, d)[2]

#=
Comparisons =#
Base.isequal(p1::P, p2::P) where {P <: AbstractPolynomial} = hash(p1) == hash(p2)
Base.:(==)(p1::AbstractPolynomial, p2::AbstractPolynomial) =
    check_same_variable(p1,p2) && (coeffs(p1) == coeffs(p2))
Base.:(==)(p::AbstractPolynomial, n::Number) = degree(p) <= 0 && constantterm(p) == n
Base.:(==)(n::Number, p::AbstractPolynomial) = p == n

function Base.isapprox(p1::AbstractPolynomial{T,X},
                       p2::AbstractPolynomial{S,Y};
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {T,X,S,Y}
    assert_same_variable(p1, p2)
    (hasnan(p1) || hasnan(p2)) && return false  # NaN poisons comparisons
    # copy over from abstractarray.jl
    Δ  = norm(p1-p2)
    if isfinite(Δ)
        return Δ <= max(atol, rtol*max(norm(p1), norm(p2)))
    else
        for i in 0:max(degree(p1), degree(p2))
            isapprox(p1[i], p2[i]; rtol=rtol, atol=atol) || return false
        end
        return true
    end
end

function Base.isapprox(p1::AbstractPolynomial{T},
                       n::S;
                       rtol::Real = (Base.rtoldefault(T, S, 0)),
                       atol::Real = 0,) where {T,S}
    return isapprox(p1, _convert(p1, [n]))
end

Base.isapprox(n::S,
    p1::AbstractPolynomial{T};
    rtol::Real = (Base.rtoldefault(T, S, 0)),
    atol::Real = 0,) where {T,S} = isapprox(p1, n, rtol = rtol, atol = atol)
