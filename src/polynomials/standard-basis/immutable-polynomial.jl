## Immutable dense / standard basis specific polynomial code
"""
    ImmutablePolynomial{T, X, N}(coeffs)

Construct an immutable (static) polynomial from its coefficients
`a₀, a₁, …, aₙ`,
lowest order first, optionally in terms of the given variable `x`
where `x` can be a character, symbol, or string.

If ``p = a_n x^n + \\ldots + a_2 x^2 + a_1 x + a_0``, we construct
this through `ImmutablePolynomial((a_0, a_1, ..., a_n))` (assuming
`a_n ≠ 0`). As well, a vector or number can be used for construction.


The usual arithmetic operators are overloaded to work with polynomials
as well as with combinations of polynomials and scalars. However,
operations involving two non-constant polynomials of different variables causes an
error. Unlike other polynomials, `setindex!` is not defined for `ImmutablePolynomials`.

As the degree of the polynomial (`+1`) is a compile-time constant,
several performance improvements are possible. For example, immutable
polynomials can take advantage of faster polynomial evaluation
provided by `evalpoly` from Julia 1.4; similar methods are also used
for addition and multiplication.

However, as the degree is included in the type, promotion between
immutable polynomials can not promote to a common type. As such, they
are precluded from use in rational functions.

!!! note
    `ImmutablePolynomial` is not axis-aware, and it treats `coeffs` simply as a
    list of coefficients with the first index always corresponding to the
    constant term.

# Examples

```jldoctest
julia> using Polynomials

julia> ImmutablePolynomial((1, 0, 3, 4))
ImmutablePolynomial(1 + 3*x^2 + 4*x^3)

julia> ImmutablePolynomial((1, 2, 3), :s)
ImmutablePolynomial(1 + 2*s + 3*s^2)

julia> one(ImmutablePolynomial)
ImmutablePolynomial(1.0)
```

!!! note
    This was modeled after
    [StaticUnivariatePolynomials](https://github.com/tkoolen/StaticUnivariatePolynomials.jl)
    by `@tkoolen`.

"""
const ImmutablePolynomial = ImmutableDensePolynomial{StandardBasis}
export ImmutablePolynomial

_typealias(::Type{P}) where {P<:ImmutablePolynomial} = "ImmutablePolynomial"

evalpoly(x, p::ImmutableDensePolynomial{B,T,X,0}) where {B<:StandardBasis,T,X} = zero(T)*zero(x)
evalpoly(x, p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,N} = EvalPoly.evalpoly(x, p.coeffs)

constantterm(p::ImmutableDensePolynomial{B,T,X,0}) where {B<:StandardBasis,T,X} = zero(T)
constantterm(p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,N} = p.coeffs[1]

scalar_add(c::S, p::ImmutableDensePolynomial{B,T,X,0}) where {B<:StandardBasis,T,X,S} =
    ImmutableDensePolynomial{B,promote_type(T,S),X,1}((c,))
function scalar_add(c::S, p::ImmutableDensePolynomial{B,T,X,1}) where {B<:StandardBasis,T,X,S}
    R = promote_type(T,S)
    ImmutableDensePolynomial{B,R,X,1}(NTuple{1,R}(p[0] + c))
end
function scalar_add(c::S, p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,S,N}
    R = promote_type(T,S)
    P = ImmutableDensePolynomial{B,R,X}
    iszero(c) && return P{N}(convert(NTuple{N,R}, p.coeffs))

    cs = _tuple_combine(+, convert(NTuple{N,R}, p.coeffs), NTuple{1,R}((c,)))
    q = P{N}(cs)

    return q
end


# return N*M
# intercept promotion call
# function Base.:*(p::ImmutableDensePolynomial{B,T,X,N},
#                  q::ImmutableDensePolynomial{B,S,X,M}) where {B<:StandardBasis,T,S,X,N,M}
#     ⊗(p,q)
# end

Base.:*(p::ImmutableDensePolynomial{B,T,X,0},
  q::ImmutableDensePolynomial{B,S,X,M}) where {B<:StandardBasis,T,S,X,M} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})
Base.:*(p::ImmutableDensePolynomial{B,T,X,N},
  q::ImmutableDensePolynomial{B,S,X,0}) where {B<:StandardBasis,T,S,X,N} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})
Base.:*(p::ImmutableDensePolynomial{B,T,X,0},
  q::ImmutableDensePolynomial{B,S,X,0}) where {B<:StandardBasis,T,S,X} = zero(ImmutableDensePolynomial{B,promote_type(T,S),X,0})
function Base.:*(p::ImmutableDensePolynomial{B,T,X,N},
                 q::ImmutableDensePolynomial{B,S,X,M}) where {B<:StandardBasis,T,S,X,N,M}
    cs = fastconv(p.coeffs, q.coeffs)
    R = eltype(cs)
    ImmutableDensePolynomial{B,R,X,N+M-1}(cs)
end


# This can be *really* slow using power_by_squaring the first time. Here we trade off a bit
# julia> p = ImmutablePolynomial(-5:5);

# julia> @time p^15;
#   0.412306 seconds (502.65 k allocations: 34.132 MiB, 6.21% gc time, 99.95% compilation time)

# julia> @time p^15;
#   0.000023 seconds (21 allocations: 5.547 KiB)

# julia> @time Base.power_by_squaring(p,15);
#  43.284660 seconds (20.41 M allocations: 1.013 GiB, 1.31% gc time, 100.00% compilation time)

# julia> @time Base.power_by_squaring(p,15);
#   0.000145 seconds (7 allocations: 6.547 KiB)
# This is not inferable, as `n` is not a compile time constant
Base.:^(p::ImmutablePolynomial, n::Integer) = immutable_power(p, n)
function immutable_power(p::ImmutablePolynomial{T,X,N}, n::Integer) where {T,X,N}
    iszero(p) && return p
    isone(N) && return ImmutablePolynomial{T,X,1}(p[0]^n)
    qs =  (PnPolynomial(p)^n).coeffs
    m = length(qs)
    N′ = n * (N-1) + 1
    z = zero(p[0])
    ImmutablePolynomial{T,X,N′}(ntuple(i -> i ≤ m ? qs[i] : z, Val(N′)))
end

#
function polynomial_composition(p::ImmutableDensePolynomial{B,T,X,N}, q::ImmutableDensePolynomial{B,S,X,M}) where {B<:StandardBasis,T,S,X,N,M}
    P = ImmutableDensePolynomial{B,promote_type(T,S), X, N*M}
    cs = evalpoly(q, p.coeffs)
    convert(P, cs)
end
function polynomial_composition(p::AbstractUnivariatePolynomial{B,T,X}, q::ImmutableDensePolynomial{B,S,Y,0}) where {B<:StandardBasis,T,S,X,Y}
    P = ImmutableDensePolynomial{B,promote_type(T,S), Y,0}
    zero(P)
end

function polynomial_composition(p::AbstractUnivariatePolynomial{B,T,X}, q::ImmutableDensePolynomial{B,S,Y,1}) where {B<:StandardBasis,T,S,X,Y}
    P = ImmutableDensePolynomial{B,promote_type(T,S), Y,1}
    P(evalpoly(constantterm(q),p))
end

# special cases of polynomial composition

# ... TBD ...


# special cases are much more performant
derivative(p::ImmutableDensePolynomial{B,T,X,0}) where {B<:StandardBasis,T,X} = p
derivative(p::ImmutableDensePolynomial{B,T,X,1}) where {B<:StandardBasis,T,X} = zero(p)

function derivative(p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,N}
    hasnan(p) && return ImmutableDensePolynomial{B,T,X,1}(zero(T)/zero(T)) # NaN{T}
    cs = ntuple(i -> i*p.coeffs[i+1], Val(N-1))
    R = eltype(cs)
    ImmutableDensePolynomial{B,R,X,N-1}(cs)
end

integrate(p::ImmutableDensePolynomial{B, T,X,0}) where {B<:StandardBasis,T,X} =
    ImmutableDensePolynomial{B,Base.promote_op(/,T,Int),X,1}((0/1,))
function integrate(p::ImmutableDensePolynomial{B,T,X,N}) where {B<:StandardBasis,T,X,N}
    N == 0 && return p # different type
    R′ = Base.promote_op(/,T,Int)
    hasnan(p) && return ImmutableDensePolynomial{B,R′,X,1}(zero(T)/zero(T)) # NaN{T}
    z = zero(first(p.coeffs))
    cs = ntuple(i -> i > 1 ? p.coeffs[i-1]/(i-1) : z/1, Val(N+1))
    R = eltype(cs)
    ImmutableDensePolynomial{B,R,X,N+1}(cs)
end
