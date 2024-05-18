using LinearAlgebra

export fromroots,
       truncate!,
       chop!,
       coeffs,
       degree,
       mapdomain,
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
```jldoctest
julia> r = [3, 2]; # (x - 3)(x - 2)

julia> fromroots(r)
Polynomial(6 - 5*x + x^2)
```
"""
function fromroots(P::Type{<:AbstractPolynomial}, rs; var::SymbolLike = :x)
    x = variable(P, var)
    p = prod(x-r for r ∈ rs; init=one(x))
    p = truncate!!(p)
    p
end
fromroots(r::AbstractVector{<:Number}; var::SymbolLike = :x) =
    fromroots(Polynomial, r, var = var)

"""
    fromroots(::AbstractMatrix{<:Number}; var=:x)
    fromroots(::Type{<:AbstractPolynomial}, ::AbstractMatrix{<:Number}; var=:x)

Construct a polynomial of the given type using the eigenvalues of the given matrix as the roots. If no type is given, defaults to `Polynomial`.

# Examples
```jldoctest
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
Vᵢⱼ βⱼ|² = argmin_β || √(W)⋅(y - V(x)β)||²`, where, `W` the diagonal
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

For fitting with a large degree, the Vandermonde matrix is exponentially
ill-conditioned. The `ArnoldiFit` type introduces an Arnoldi orthogonalization
that fixes this problem.


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
fit′(P::Type{<:AbstractPolynomial}, x, y, args...;kwargs...) = throw(ArgumentError("x and y do not produce abstract vectors"))
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
             deg = length(x) - 1;
             weights = nothing,
             var = :x,) where {T}
    x = mapdomain(P, x)
    vand = vander(P, x, deg)
    if !isnothing(weights)
        coeffs = _wlstsq(vand, y, weights)
    else
        coeffs = vand \ y
    end
    R = float(T)
    if isa(deg, Integer)
        return P(R.(coeffs), var)
    else
        cs = zeros(T, 1 + maximum(deg))
        for (i,aᵢ) ∈ zip(deg, coeffs)
            cs[1 + i] = aᵢ
        end
        return P(cs, var)
    end


end


# Weighted linear least squares
_wlstsq(vand, y, W::Number) = _wlstsq(vand, y, fill!(similar(y), W))
function _wlstsq(vand, y, w::AbstractVector)
    W = Diagonal(sqrt.(w))
    (W * vand) \ (W * y)
end
_wlstsq(vand, y, W::AbstractMatrix) = (vand' * W * vand) \ (vand' * W * y)

"""
    roots(::AbstractPolynomial; kwargs...)

Returns the roots, or zeros, of the given polynomial.

For non-factored, standard basis polynomials the roots are calculated via the
eigenvalues of the companion matrix. The `kwargs` are passed to the
`LinearAlgebra.eigvals` call.

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

Returns the indefinite integral of the polynomial with constant `C` when
expressed in the standard basis.
"""
function integrate(p::P, C) where {P <: AbstractPolynomial}
    ∫p = integrate(p)
    isnan(C) && return ⟒(P){eltype(∫p+C), indeterminate(∫p)}([C])
    ∫p + (C - constantterm(∫p))
end

"""
    integrate(::AbstractPolynomial, a, b)

Compute the definite integral of the given polynomial from `a` to `b`. Will
throw an error if either `a` or `b` are out of the polynomial's domain.
"""
function integrate(p::AbstractPolynomial, a, b)
    P = integrate(p)
    return P(b) - P(a)
end

"""
    derivative(::AbstractPolynomial, order::Int = 1)

Returns a polynomial that is the `order`th derivative of the given polynomial.
`order` must be non-negative.
"""
derivative(::AbstractPolynomial, ::Int)


"""
    critical_points(p::AbstractPolynomial{<:Real}, I=domain(p); endpoints::Bool=true)

Return the critical points of `p` (real zeros of the derivative) within `I` in sorted order.

* `p`: a polynomial

* `I`: a specification of a closed or infinite domain, defaulting to
  `Polynomials.domain(p)`. When specified, the values of `extrema(I)` are used
  with closed endpoints when finite.

* `endpoints::Bool`: if `true`, return the endpoints of `I` along with the critical points


Can be used in conjunction with `findmax`, `findmin`, `argmax`, `argmin`, `extrema`, etc.

## Example
```julia
x = variable()
p = x^2 - 2
cps = Polynomials.critical_points(p)
findmin(p, cps)  # (-2.0, 2.0)
argmin(p, cps)   #  0.0
extrema(p, cps)  # (-2.0, Inf)
cps = Polynomials.critical_points(p, (0, 2))
extrema(p, cps)  # (-2.0, 2.0)
```

!!! note
    There is a *big* difference between `minimum(p)` and `minimum(p, cps)`. The former takes the viewpoint that a polynomial `p` is a certain type of vector of its coefficients; returning the smallest coefficient. The latter uses `p` as a callable object, returning the smallest of the values `p.(cps)`.
"""
function critical_points(p::AbstractPolynomial{T}, I = domain(p);
                         endpoints::Bool=true) where {T <: Real}

    I′ = Interval(I)
    l, r = extrema(I′)

    q = Polynomials.ngcd(derivative(p), derivative(p,2)).v
    pts = sort(real.(filter(isreal, roots(q))))
    pts = filter(in(I′), pts)

    !endpoints && return pts

    l !== first(pts) && pushfirst!(pts, l)
    r != last(pts) && push!(pts, r)
    pts
end


## --------------------------------------------------

"""
    truncate!(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

In-place version of [`truncate`](@ref)
"""
truncate!(p::AbstractPolynomial; kwargs...) = _truncate!(p; kwargs...)

function _truncate!(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {T}
    _truncate!(p.coeffs, rtol=rtol, atol=atol)
    chop!(p, rtol = rtol, atol = atol)
end

## _truncate! underlying storage type
function _truncate!(ps::Vector{T};
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

function _truncate!(ps::Dict{S,T};
                   rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {S,T}

    isempty(ps) && return nothing
    max_coeff = norm(values(ps), Inf)
    thresh = max_coeff * rtol + atol

    for (k,val) in  ps
        if abs(val) <= thresh
            pop!(ps,k)
        end
    end
    nothing
end

_truncate!(ps::NTuple; kwargs...) = throw(ArgumentError("`truncate!` not defined for tuples."))

# _truncate(ps::NTuple{0}; kwargs...) = ps
# function _truncate(ps::NTuple{N,T};
#                    rtol::Real = Base.rtoldefault(real(T)),
#                    atol::Real = 0,) where {N,T}


#     thresh = norm(ps, Inf) * rtol + atol
#     return NTuple{N,T}(abs(pᵢ) <= thresh ? zero(T) : pᵢ for pᵢ ∈ values(ps))
# end


"""
    truncate(::AbstractPolynomial{T};
        rtol::Real = Base.rtoldefault(real(T)), atol::Real = 0)

Rounds off coefficients close to zero, as determined by `rtol` and `atol`, and then chops any leading zeros. Returns a new polynomial.
"""
function Base.truncate(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
    atol::Real = 0,) where {T}
    _truncate!(deepcopy(p), rtol = rtol, atol = atol)
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
    isempty(ps) && return ps
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

_chop(ps::NTuple{0}; kwargs...) = ps
function _chop(ps::NTuple{N};
               rtol = nothing,
               atol = nothing) where {N}

    T = real(eltype(ps))
    rtol = something(rtol, Base.rtoldefault(T))
    atol = something(atol, zero(T))

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

Removes any leading coefficients that are approximately 0 (using `rtol` and `atol` with `norm(p)`). Returns a polynomial whose degree will guaranteed to be equal to or less than the given polynomial's.
"""
function Base.chop(p::AbstractPolynomial{T};
    rtol::Real = Base.rtoldefault(real(T)),
    atol::Real = 0,) where {T}
    chop!(deepcopy(p), rtol = rtol, atol = atol)
end


# for generic usage, as immutable types are not mutable
chop!!(p::AbstractPolynomial; kwargs...) = (p = chop!(p; kwargs...); p)
truncate!!(p::AbstractPolynomial; kwargs...) = _truncate!(p; kwargs...)

## --------------------------------------------------

"""
    check_same_variable(p::AbstractPolynomial, q::AbstractPolynomial)

Check if either `p` or `q` is constant or if `p` and `q` share the same variable
"""
check_same_variable(p::AbstractPolynomial, q::AbstractPolynomial) =
    (isconstant(p) || isconstant(q)) || indeterminate(p) ==  indeterminate(q)

function assert_same_variable(p::AbstractPolynomial, q::AbstractPolynomial)
    check_same_variable(p,q) || throw(ArgumentError("Non-constant polynomials have different indeterminates"))
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
function LinearAlgebra.norm(q::AbstractPolynomial{T,X}, p::Real = 2) where {T,X}
    iszero(q) && return zero(real(T))^(1/p)
    return norm(values(q), p)
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
function Base.convert(P::Type{<:AbstractPolynomial}, q::AbstractPolynomial)
    X = indeterminate(P,q)
    x = variable(P, X)
    q(x)
end
function Base.convert(::Type{T}, p::AbstractPolynomial{T,X}) where {T <: Number,X}
    isconstant(p) && return T(constantterm(p))
    throw(ArgumentError("Can't convert a nonconstant polynomial to type $T"))
end



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
Base.size(p::AbstractPolynomial) = (length(p),)
Base.size(p::AbstractPolynomial, i::Integer) =  i <= 1 ? size(p)[i] : 1
Base.eltype(p::AbstractPolynomial{T}) where {T} = T
# in  analogy  with  polynomial as a Vector{T} with different operations defined.
Base.eltype(::Type{<:AbstractPolynomial}) = Float64
Base.eltype(::Type{<:AbstractPolynomial{T}}) where {T} = T
_eltype(::Type{<:AbstractPolynomial}) = nothing
_eltype(::Type{<:AbstractPolynomial{T}}) where {T} = T
function _eltype(P::Type{<:AbstractPolynomial}, p::AbstractPolynomial)
    T′ = _eltype(P)
    T = isnothing(T′) ? eltype(p) : T′
    T
end

"""
    copy_with_eltype(::Type{T}, [::Val{X}], p::AbstractPolynomial)

Copy polynomial `p` changing the underlying element type and optionally the symbol.
"""
copy_with_eltype(::Type{T}, ::Val{X}, p::P) where {T, X, S, Y, P <:AbstractPolynomial{S,Y}} =
    ⟒(P){T, Symbol(X)}(p.coeffs)
copy_with_eltype(::Type{T}, p::P) where {T, S, Y, P <:AbstractPolynomial{S,Y}} =
    copy_with_eltype(T, Val(Y), p)
# easier to type if performance isn't an issue, but could be dropped
#copy_with_eltype(::Type{T}, X, p::P) where {T, S, Y, P<:AbstractPolynomial{S, Y}} =
#    copy_with_eltype(Val(T), Val(X), p)
#copy_with_eltype(::Type{T}, p::P) where {T, S, X, P<:AbstractPolynomial{S,X}} =
#    copy_with_eltype(Val(T), Val(X), p)

"""
    iszero(p::AbstractPolynomial)

Is this a ``0`` polynomial.

For most types, the ``0`` polynomial is one with no coefficients (coefficient vector `T[]`),
though some types have the possibility of trailing zeros. The degree of a zero polynomial is conventionally ``-1``, though this is not the convention for Laurent polynomials.
"""
Base.iszero(p::AbstractPolynomial) = all(iszero, values(p))::Bool


# See discussions in https://github.com/JuliaMath/Polynomials.jl/issues/258
"""
    all(pred, poly::AbstractPolynomial)

Test whether all coefficients of an `AbstractPolynomial` satisfy predicate `pred`.

You can implement `isreal`, etc., to a `Polynomial` by using `all`.
"""
Base.all(pred, p::AbstractPolynomial{T, X}) where {T,X} = all(pred, values(p))
"""
    any(pred, poly::AbstractPolynomial)

Test whether any coefficient of an `AbstractPolynomial` satisfies predicate `pred`.
"""
Base.any(pred, p::AbstractPolynomial{T,X}) where {T, X} = any(pred, values(p))




"""
    map(fn, p::AbstractPolynomial, args...)

Transform coefficients of `p` by applying a function (or other callables) `fn` to each of them.

You can implement `real`, etc., to a `Polynomial` by using `map`. The type of `p` may narrow using this function.
"""
function Base.map(fn, p::P, args...)  where {P<:AbstractPolynomial}
    xs = map(fn, p.coeffs, args...)
    R = eltype(xs)
    X = indeterminate(p)
    return ⟒(P){R,X}(xs)
end


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
hasnan(p::AbstractPolynomial) = any(hasnan, p)::Bool
hasnan(p::AbstractArray) = any(hasnan.(p))
hasnan(x) = isnan(x)

"""
    isconstant(::AbstractPolynomial)

Is the polynomial  `p` a constant.
"""
isconstant(p::AbstractPolynomial) = degree(p) <= 0 && firstindex(p) == 0

# XXX docstring isn't quite right!
# coeffs returns
# * p.coeffs for Laurent types (caller needs to be aware that offset is possible)
# * a container with (a₀, a₁, …, aₙ) which for some types may have trailing zeros (ImmutableDense)
"""
    coeffs(::AbstractPolynomial)
    coeffs(::AbstractDenseUnivariatePolynomial)
    coeffs(::AbstractLaurentUnivariatePolynomial)

For a dense, univariate polynomial return the coefficients ``(a_0, a_1, \\dots, a_n)``
as an interable. This may be a vector or tuple, and may alias the
polynomials coefficients.

For a Laurent type polynomial (e.g. `LaurentPolynomial`, `SparsePolynomial`) return the coefficients ``(a_i, a_{i+1}, \\dots, a_j)`` where
``i`` is found from `firstindex(p)` and ``j`` from `lastindex(p)`.

For `LaurentPolynomial` and `SparsePolynomial`, the `pairs` iterator is more generically useful, as it iterates over ``(i, p_i)`` possibly skipping the terms where ``p_i = 0``.

Defaults to `p.coeffs`.
"""
coeffs(p::AbstractPolynomial) = p.coeffs



# specialize this to p[0] when basis vector is 1
"""
    constantterm(p::AbstractPolynomial)

return `p(0)`, the constant term in the standard basis
"""
constantterm(p::AbstractPolynomial{T}) where {T} = p(zero(T))

"""
    degree(::AbstractPolynomial)

Return the degree of the polynomial, i.e. the highest exponent in the polynomial that
has a nonzero coefficient.

For standard basis polynomials the degree of the zero polynomial is defined to be ``-1``.
For Laurent type polynomials, this is `0`, or `lastindex(p)`. The `firstindex` method gives the smallest power
of the indeterminate for the polynomial.
The default method assumes the basis polynomials, `βₖ`, have degree `k`.

"""
degree(p::AbstractPolynomial) = iszero(coeffs(p)) ? -1 : length(coeffs(p)) - 1 + min(0, minimumexponent(p))



"""
    Polynomials.domain(::Type{<:AbstractPolynomial})

Returns the domain of the polynomial.
"""
domain(::Type{<:AbstractPolynomial})
domain(::P) where {P <: AbstractPolynomial} = domain(P)

"""
    mapdomain(::Type{<:AbstractPolynomial}, x::AbstractArray)
    mapdomain(::AbstractPolynomial, x::AbstractArray)

Given values of x that are assumed to be unbounded (-∞, ∞), return values rescaled to the domain of the given polynomial.

# Examples
```jldoctest
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
# minimumexponent(p) returns min(0, minimum(keys(p)))
# For most polynomial types, this is statically known to be zero
# For polynomials that support negative indices, minimumexponent(typeof(p))
# should return typemin(Int)
minimumexponent(p::AbstractPolynomial) = minimumexponent(typeof(p))
minimumexponent(::Type{<:AbstractPolynomial}) = 0
# firstindex, lastindex correspond to the range of which basis vectors are being represented by the coefficients
"""
    firstindex(p::AbstractPolynomial)

The index of the smallest basis element, ``\\beta_i``,  represented by the coefficients. This is ``0`` for
a zero polynomial.
"""
Base.firstindex(p::AbstractPolynomial) = 0  # XXX() is a better default
"""
    lastindex(p::AbstractPolynomial)

The index of the largest basis element, ``\\beta_i``,  represented by the coefficients.
May be ``-1`` or ``0`` for the zero polynomial, depending on the storage type.
"""
Base.lastindex(p::AbstractPolynomial) = length(p) - 1 + firstindex(p) # not degree, which accounts for any trailing zeros
"""
    eachindex(p::AbstractPolynomial)

Iterator over all indices of the represented basis elements
"""
Base.eachindex(p::AbstractPolynomial) = firstindex(p):lastindex(p)
Base.broadcastable(p::AbstractPolynomial) = Ref(p)
degreerange(p::AbstractPolynomial) = firstindex(p):lastindex(p)

# getindex
function Base.getindex(p::AbstractPolynomial{T}, idx::Int) where {T}
    m,M = firstindex(p), lastindex(p)
    m > M && return zero(T)
    idx < m && throw(BoundsError(p, idx))
    idx > M && return zero(T)
    p.coeffs[idx - m + 1]
end
#XXXBase.getindex(p::AbstractPolynomial, idx::Number) = getindex(p, convert(Int, idx))
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
Base.setindex!(p::AbstractPolynomial, values, indices) =
    [setindex!(p, v, i) for (v, i) in tuple.(values, indices)]
Base.setindex!(p::AbstractPolynomial, values, ::Colon) =
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
function _iterate(p, state)
    firstindex(p) <= state <= lastindex(p) || return nothing
    return p[state], state+1
end
Base.iterate(p::AbstractPolynomial, state = firstindex(p)) = _iterate(p, state)

# pairs map i -> aᵢ *possibly* skipping over ai == 0
# cf. abstractdict.jl
struct PolynomialKeys{P} <: AbstractSet{Int}
    p::P
end
struct PolynomialValues{P, T}
    p::P

    PolynomialValues{P}(p::P) where {P} = new{P, eltype(p)}(p)
    PolynomialValues(p::P) where {P} = new{P, eltype(p)}(p)
end
"""
    keys(p::AbstractPolynomial)

Iterator over ``i``s for each basis element, ``\\beta_i``, represented by the coefficients.
"""
Base.keys(p::AbstractPolynomial) =  PolynomialKeys(p)
"""
    values(p::AbstractPolynomial)

Iterator over ``p_i``s for each basis element, ``\\beta_i``, represented by the coefficients.
"""
Base.values(p::AbstractPolynomial) =  PolynomialValues(p)
"""
    pairs(p::AbstractPolynomial)

Iterator over ``(i, p_i)`` for each basis element, ``\\beta_i``, represented by the coefficients.
"""
Base.pairs(p::AbstractPolynomial)
Base.length(p::PolynomialValues) = length(p.p.coeffs)
Base.eltype(p::PolynomialValues{<:Any,T}) where {T} = T
Base.length(p::PolynomialKeys) = length(p.p.coeffs)
Base.size(p::Union{PolynomialValues, PolynomialKeys}) = (length(p),)
function Base.iterate(v::PolynomialKeys, state = firstindex(v.p))
    firstindex(v.p) <= state <= lastindex(v.p) || return nothing
    return state, state+1
end

Base.iterate(v::PolynomialValues, state = firstindex(v.p)) = _iterate(v.p, state)


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
    isnothing(y) && return nothing
    kv, s = y
    return (kv[2]*basis(v.p, kv[1]), s)
end
Base.length(v::Monomials) = length(keys(v.p))


#=
identity =#
Base.copy(p::P) where {P <: AbstractPolynomial} = map(identity, p)
Base.hash(p::AbstractPolynomial{T,X}, h::UInt) where {T,X} = hash(indeterminate(p), hash(p.coeffs, hash(X,h)))

# get symbol of polynomial. (e.g. `:x` from 1x^2 + 2x^3...
_indeterminate(::Type{P}) where {P <: AbstractPolynomial} = nothing
_indeterminate(::Type{P}) where {T, X, P <: AbstractPolynomial{T,X}} = X
indeterminate(::Type{P}) where {P <: AbstractPolynomial} = something(_indeterminate(P), :x)
indeterminate(p::P) where {P <: AbstractPolynomial} = _indeterminate(P)

function indeterminate(PP::Type{P}, p::AbstractPolynomial{T,Y}) where {P <: AbstractPolynomial, T,Y}
    X = _indeterminate(PP)
    isnothing(X) && return Y
    isconstant(p) && return X
    assert_same_variable(X,Y)
    return X
    #X = isnothing(_indeterminate(PP)) ? indeterminate(p) :  _indeterminate(PP)
end
indeterminate(PP::Type{P}, x::Symbol) where {P <: AbstractPolynomial} = something(_indeterminate(PP), x)

#=
zero, one, variable, basis =#



"""
    zero(::Type{<:AbstractPolynomial})
    zero(::AbstractPolynomial)

Returns a representation of 0 as the given polynomial.
"""
function Base.zero(::Type{P}) where {P<:AbstractPolynomial}
    T,X = eltype(P), indeterminate(P)
    ⟒(P){T,X}(T[])
end
Base.zero(::Type{P}, var::SymbolLike) where {P <: AbstractPolynomial} = zero(⟒(P){eltype(P),Symbol(var)}) #default 0⋅b₀
Base.zero(p::P, var=indeterminate(p)) where {P <: AbstractPolynomial} = zero(P, var)
"""
    one(::Type{<:AbstractPolynomial})
    one(::AbstractPolynomial)

Returns a representation of 1 as the given polynomial.
"""
Base.one(::Type{P}) where {P <: AbstractPolynomial} =  one(⟒(P){eltype(P), indeterminate(P)})
Base.one(::Type{P}, var::SymbolLike) where {P <: AbstractPolynomial} = one(⟒(P){eltype(P), Symbol(var)})
Base.one(p::P, var=indeterminate(p)) where {P <: AbstractPolynomial} = one(P, var)
# each polynomial type implements:
# Base.one(::Type{P}) where {T,X,P<:AbstractPolynomial{T,X}} = throw(ArgumentError("No default method defined")) # no default method

Base.oneunit(::Type{P}, args...) where {P <: AbstractPolynomial} = one(P, args...)
Base.oneunit(p::P, args...) where {P <: AbstractPolynomial} = one(p, args...)


"""
    variable(var=:x)
    variable(::Type{<:AbstractPolynomial}, var=:x)
    variable(p::AbstractPolynomial, var=indeterminate(p))

Return the monomial `x` in the indicated polynomial basis.  If no type is give, will default to [`Polynomial`](@ref). Equivalent to `P(var)`.

# Examples
```jldoctest
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
variable(::Type{P}) where {P <: AbstractPolynomial} = variable(⟒(P){eltype(P), indeterminate(P)})
variable(::Type{P}, var::SymbolLike) where {P <: AbstractPolynomial} = variable(⟒(P){eltype(P),Symbol(var)})
variable(p::AbstractPolynomial, var = indeterminate(p)) = variable(⟒(p){eltype(p), Symbol(var)})
variable(var::SymbolLike = :x) = variable(Polynomial{Int,Symbol(var)})
# Each polynomial type implements:
# variable(::Type{P}) where {T,X,P <: AbstractPolynomial{T,X}} = throw(ArgumentError("No default method defined")) # no default

# Exported in #470. Exporting was a mistake!
# Now must be used through qualification: `Polynomials.@variable ...`
#@variable x
#@variable x::Polynomial
#@variable x::Polynomial{t]
macro variable(x)
    q = Expr(:block)
    if isa(x, Expr) && x.head == :(::)
        x, P = x.args
        push!(q.args, Expr(:(=), esc(x),
                           Expr(:call, :variable, P, Expr(:quote, x))))
    else
        push!(q.args, Expr(:(=), esc(x),
                           Expr(:call, :variable, Expr(:quote, x))))
    end
    push!(q.args, esc(x))

    q
end



# basis
# var is a positional argument, not a keyword; can't deprecate so we do `_var; var=_var`
# return the kth basis polynomial for the given polynomial type, e.g. x^k for Polynomial{T}
"""
    basis(p::P, i::Int)
    basis(::Type{<:AbstractPolynomial}, i::Int, var=:x)

Return ith basis element for a given polynomial type, optionally with a specified variable.
"""
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
composition
cf. https://github.com/JuliaMath/Polynomials.jl/issues/511 for a paper with implementations

=#
"""
    polynomial_composition(p, q)

Evaluate `p(q)`, possibly exploiting a faster evaluation scheme, defaulting to `evalpoly`.
"""
function polynomial_composition(p::AbstractPolynomial, q::AbstractPolynomial)
    evalpoly(q, p)
end

#=
arithmetic =#
const Scalar = Union{Number, Matrix}

# scalar operations
# scalar_add utilized polynomial addition. May be more performant to provide new method
scalar_add(c::S, p::AbstractPolynomial) where {S} = p + c * one(p)

# Scalar multiplication; no assumption of commutivity
scalar_mult(p::P, c::S) where {S, T, X, P<:AbstractPolynomial{T,X}} = map(Base.Fix2(*,c), p)
scalar_mult(c::S, p::P) where {S, T, X, P<:AbstractPolynomial{T,X}} = map(Base.Fix1(*,c), p)
scalar_mult(p1::AbstractPolynomial, p2::AbstractPolynomial) =
    throw(ArgumentError("scalar_mult(::$(typeof(p1)), ::$(typeof(p2))) is not defined.")) # avoid ambiguity, issue #435

# scalar div (faster than Base.Fix2(/,c) )
scalar_div(p::P, c::S) where {S, T, X, P<:AbstractPolynomial{T, X}} = map(Base.Fix2(*, one(T)/c), p)



Base.:-(p::P) where {P <: AbstractPolynomial} = map(-, p)

Base.:*(p::AbstractPolynomial, c::Scalar) = scalar_mult(p, c)
Base.:*(c::Scalar, p::AbstractPolynomial) = scalar_mult(c, p)
Base.:*(c::T, p::P) where {T, X, P <: AbstractPolynomial{T,X}} = scalar_mult(c, p)
Base.:*(p::P, c::T) where {T, X, P <: AbstractPolynomial{T,X}} = scalar_mult(p, c)

# implicitly identify c::Scalar with a constant polynomials
Base.:+(c::Scalar, p::AbstractPolynomial) = scalar_add(c, p)
Base.:-(p::AbstractPolynomial, c::Scalar) = scalar_add(-c, p)
Base.:-(c::Scalar, p::AbstractPolynomial) = scalar_add(c, -p) # extra copy! eww

## addition
## Fall back addition is possible as vector addition with padding by 0s
## Subtypes will likely want to implement both:
## +(p::P,c::Scalar) and +(p::P, q::Q) where {T,S,X,P<:SubtypePolynomial{T,X},Q<:SubtypePolynomial{S,X}}
## though the default for poly+poly isn't terrible

Base.:+(p::AbstractPolynomial) = p

# polynomial + scalar; implicit identification of c with c*one(p)
# what are these here for?
Base.:+(p::P, c::T) where {T,X, P<:AbstractPolynomial{T,X}} = scalar_add(c, p)
Base.:+(p::P, c::S) where {T,X, P<:AbstractPolynomial{T,X}, S} = scalar_add(c,p)

## polynomial + polynomial when different types
## - each polynomial container type implements PB{B,T,X} + PB{B,S,X}
## - this handles case X ≠ Y unless constant
## - when PB₁ ≠ PB₂ we promote both polynomials
function Base.:+(p::P, q::Q) where {T,X,P <: AbstractPolynomial{T,X}, S,Y,Q <: AbstractPolynomial{S,Y}}
    isconstant(p) && return constantterm(p) + q
    isconstant(q) && return p + constantterm(q)
    assert_same_variable(X,Y) # should error
end

function Base.:+(p::P, q::Q) where {T,X,P <: AbstractPolynomial{T,X}, S,Q <: AbstractPolynomial{S,X}}
    sum(promote(p,q))
end


function Base.:-(p::P, q::Q) where {T,X,P <: AbstractPolynomial{T,X}, S,Y,Q <: AbstractPolynomial{S,Y}}
    isconstant(p) && return constantterm(p) + q
    isconstant(q) && return p + constantterm(q)
    assert_same_variable(X,Y) # should error
end

function Base.:-(p::P, q::Q) where {T,X,P <: AbstractPolynomial{T,X}, S,Q <: AbstractPolynomial{S,X}}
    -(promote(p,q)...)
end

# the case p::P{B,T,X} - q::P{B,S,X} should be done in container types
Base.:-(p::P, q::P) where {T,X,P <: AbstractPolynomial{T,X}} = p + (-1*q)


## -- multiplication
## Polynomial p*q
## Polynomial multiplication formula depend on the particular basis used.
## The subtype must implement *(::PT{T,X,[N]}, ::PT{S,X,[M]})
function Base.:*(p1::P, p2::Q) where {T,X,P <: AbstractPolynomial{T,X},
                                      S,Y,Q <: AbstractPolynomial{S,Y}}
    isconstant(p1) && return constantterm(p1) * p2
    isconstant(p2) && return p1 * constantterm(p2)
    assert_same_variable(X, Y) # should error
end

function Base.:*(p1::P, p2::Q) where {T,X,P <: AbstractPolynomial{T,X},
                                      S,  Q <: AbstractPolynomial{S,X}}
    prod(promote(p1, p2))
end


# scalar div
Base.:/(p::AbstractPolynomial, c) = scalar_div(p, c)

Base.:^(p::AbstractPolynomial, n::Integer) = Base.power_by_squaring(p, n)


## -----

# <<< could be moved to SpecialPolynomials
# Works when p,q of same type.
# For Immutable, must remove N,M bit;
# for others, can widen to Type{T,X}, Type{S,X} to avoid a promotion
# function Base.:+(p::P, q::P) where {T,X,P<:AbstractPolynomial{T,X}}
#     cs = degree(p) >= degree(q)  ? ⊕(P, p.coeffs, q.coeffs) : ⊕(P, q.coeffs, p.coeffs)
#     return P(cs)
# end



# Th ⊕ method below is used in Special Polynomials, but not here, as it was removed for
# similar methods in the polynomial-basetypes
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
@generated function ⊕(P::Type{<:AbstractPolynomial}, p1::NTuple{N}, p2::NTuple{M}) where {N,M}
    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(p1[$i] + p2[$i])
    end
    for i in M+1:N
        exprs[i] =:(p1[$i])
    end

    return quote
        Base.@_inline_meta
        #Base.@inline
        tuple($(exprs...))
    end

end

# addition when a dictionary is used for storage
function ⊕(P::Type{<:AbstractPolynomial}, p1::Dict{Int,T}, p2::Dict{Int,S}) where {T,S}
    R = promote_type(T,S)
    p = Dict{Int, R}()


    # this allocates in the union
    #for i in union(eachindex(p1), eachindex(p2))
    #    p[i] = p1[i] + p2[i]
    #end

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
## >>> moved to SpecialPolynomials (or rewrite that to use new container types)

## -----

function Base.divrem(num::P, den::O) where {P <: AbstractPolynomial,O <: AbstractPolynomial}
    n, d = promote(num, den)
    return divrem(n, d)
end

"""
    gcd(a::AbstractPolynomial, b::AbstractPolynomial; atol::Real=0, rtol::Real=Base.rtoldefault)

Find the greatest common denominator of two polynomials recursively using
[Euclid's algorithm](http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm).

# Examples

```jldoctest
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
#uvw(p1::P, p2::P; kwargs...) where {P <:AbstractPolynomial} = throw(ArgumentError("uvw not defined"))

"""
    div(::AbstractPolynomial, ::AbstractPolynomial)
"""
Base.div(n::AbstractPolynomial, d::AbstractPolynomial) = divrem(n, d)[1]

"""
    rem(::AbstractPolynomial, ::AbstractPolynomial)
"""
Base.rem(n::AbstractPolynomial, d::AbstractPolynomial) = divrem(n, d)[2]

#=
Comparisons
=#

Base.isequal(p1::P, p2::P) where {P <: AbstractPolynomial} = hash(p1) == hash(p2)
function Base.:(==)(p1::P, p2::P) where {P <: AbstractPolynomial}
    iszero(p1) && iszero(p2) && return true
    eachindex(p1) == eachindex(p2) || return false
    # coeffs(p1) == coeffs(p2), but non-allocating
    p1val = (p1[i] for i in eachindex(p1))
    p2val = (p2[i] for i in eachindex(p2))
    all(((a,b),) -> a == b, zip(p1val, p2val))
end
function Base.:(==)(p1::AbstractPolynomial, p2::AbstractPolynomial)
    if isconstant(p1)
        isconstant(p2) && return constantterm(p1) == constantterm(p2)
        return false
    elseif isconstant(p2)
        return false # p1 is not constant
    end
    check_same_variable(p1, p2) || return false
    ==(promote(p1,p2)...)
end
Base.:(==)(p::AbstractPolynomial, n::Scalar) = isconstant(p) && constantterm(p) == n
Base.:(==)(n::Scalar, p::AbstractPolynomial) = p == n


function Base.isapprox(p1::AbstractPolynomial, p2::AbstractPolynomial; kwargs...)
    if isconstant(p1)
        isconstant(p2) && return constantterm(p1) == constantterm(p2)
        return false
    elseif isconstant(p2)
        return false
    end
    assert_same_variable(p1, p2) || return false
    isapprox(promote(p1, p2)...; kwargs...)
end

function Base.isapprox(p1::AbstractPolynomial{T,X},
                       p2::AbstractPolynomial{S,X};
                       rtol::Real = (Base.rtoldefault(T,S,0)),
                       atol::Real = 0,) where {T,S,X}
    (hasnan(p1) || hasnan(p2)) && return false  # NaN poisons comparisons
    # copy over from abstractarray.jl
    Δ  = norm(p1-p2)
    if isfinite(Δ)
        return Δ <= max(atol, rtol*max(norm(p1), norm(p2)))
    else
        for i in minimum(firstindex, (p1,p2)):maximum(degree, (p1,p2))
            isapprox(p1[i], p2[i]; rtol=rtol, atol=atol) || return false
        end
        return true
    end
end

Base.isapprox(p1::AbstractPolynomial{T}, n::S;kwargs...) where {S,T} = isapprox(p1, n*one(p1))
Base.isapprox(n::S,p1::AbstractPolynomial{T}; kwargs...) where {S,T} = isapprox(p1, n; kwargs...)

Base.isapprox(::AbstractPolynomial{T}, ::Missing, args...; kwargs...) where T = missing
Base.isapprox(::Missing, ::AbstractPolynomial{T}, args...; kwargs...) where T = missing

function LinearAlgebra.dot(::AbstractPolynomial, ::AbstractPolynomial, args...; kwargs...)
    throw(ArgumentError("No generic `dot` method is defined for polynomials."))
end
