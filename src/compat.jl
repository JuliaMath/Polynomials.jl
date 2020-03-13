## The plan: keep these to ensure underlying changes are not disruptive
## then deprecate these
## then release v1.0

poly(r, var = :x) = fromroots(Polynomial, r; var = var)
polyval(p::AbstractPolynomial, x::Number) = p(x)
polyval(p::AbstractPolynomial, x) = p.(x)

function Base.getproperty(p::AbstractPolynomial, nm::Symbol)
    if nm == :a
        #Base.depwarn("AbstracPolynomial.a is deprecated, use AbstracPolynomial.coeffs or coeffs(AbstractPolynomial) instead.",
        #    Symbol("Base.getproperty"),
        #)
        return getfield(p, :coeffs)
    end
    return getfield(p, nm)
end

polyint(p::AbstractPolynomial, C = 0) = integrate(p, C)
polyint(p::AbstractPolynomial, a, b) = integrate(p, a, b)
polyder(p::AbstractPolynomial, ord = 1) = derivative(p, ord)
polyfit(x, y, n = length(x) - 1, sym=:x) = fit(Polynomial, x, y, n; var = sym)
polyfit(x, y, sym::Symbol) = fit(Polynomial, x, y, var = sym)

padeval(PQ::Pade, x::Number) = PQ(x)
padeval(PQ::Pade, x) = PQ.(x)

export poly, polyval, polyint, polyder, polyfit, padeval


## @deprecate poly(r, var = :x) fromroots(Polynomial, r; var = var)
## @deprecate polyval(p::AbstractPolynomial, x::Number) p(x)
## @deprecate polyval(p::AbstractPolynomial, x) p.(x)

## function Base.getproperty(p::AbstractPolynomial, nm::Symbol)
##     if nm == :a
##         Base.depwarn("AbstracPolynomial.a is deprecated, use AbstracPolynomial.coeffs or coeffs(AbstractPolynomial) instead.",
##             Symbol("Base.getproperty"),
##         )
##         return getfield(p, :coeffs)
##     end
##     return getfield(p, nm)
## end

## @deprecate polyint(p::AbstractPolynomial, C = 0) integrate(p, C)
## @deprecate polyint(p::AbstractPolynomial, a, b) integrate(p, a, b)
## @deprecate polyder(p::AbstractPolynomial, ord = 1) derivative(p, ord)
## @deprecate polyfit(x, y, n = length(x) - 1) fit(Polynomial, x, y; deg = n)
## @deprecate polyfit(x, y, sym::Symbol) fit(Polynomial, x, y; var = sym)

## @deprecate padeval(PQ::Pade, x::Number) PQ(x)
## @deprecate padeval(PQ::Pade, x) PQ.(x)
