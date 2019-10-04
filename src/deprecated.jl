@deprecate Poly Polynomial
@deprecate Poly(n::Number, var = :x) Polynomial(n, var)
@deprecate Poly(a::AbstractVector, var::SymbolLike = :x) Polynomial(a, var)
# @deprecate Poly{T}(x::AbstractVector{S}, var::SymbolLike = :x) where {T,S} Polynomial{T}(x, var)

@deprecate poly(r, var = :x) fromroots(Polynomial, r; var=var)
@deprecate polyval(p::Polynomial, x) p.(x)
@deprecate polyint(p::Polynomial, k = 0) integral(p, k)
@deprecate polyint(p::Polynomial, a, b) integrate(p, a, b)
@deprecate polyder(p::Polynomial, ord = 1) derivative(p, ord)
@deprecate polyfit(x, y, n = length(x) - 1) fit(Polynomial, x, y; deg = n)
@deprecate polyfit(x, y, sym::Symbol) fit(Polynomial, x, y; var=sym)

@deprecate padeval(PQ::Pade, x) PQ(x)
