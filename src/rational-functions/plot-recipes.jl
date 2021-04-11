## Plot recipe
## attempt to work around asymptotes
function rational_function_trim(pq, a=nothing,b=nothing)

    p,q = normal_form(convert(RationalFunction,pq), method=:numerical)
    p′,q′ = normal_form(derivative(p//q))

    λs = Polynomials.Multroot.multroot(q).values
    λs = isempty(λs) ? λs : filter(isreal, λs)

    cps = Polynomials.Multroot.multroot(p′).values
    cps = isempty(cps) ? cps : filter(isreal, cps)
    
    if a==nothing && b==nothing
        u= poly_interval(p)
        v= poly_interval(q)
        a,b = min(first(u), first(v)), max(last(u), last(v))

        if !isempty(λs)
            a,b = min(a, real(minimum(λs))), max(b, real(maximum(λs)))
        end
        if !isempty(cps)
            a,b = min(a, real(minimum(cps))), max(b, real(maximum(cps)))
        end
        a *= (a > 0 ? 1/1.5 : 1.25)
        b *= (b < 0 ? 1/1.5 : 1.25)
    end

    n = 501
    xs = range(a,stop=b, length=n)
    ys = pq.(xs)

    

    M = max(20, 3*maximum(abs, pq.(cps)))
    ys′ = [abs(y) > M ? NaN : y for y ∈ ys]

    xs, ys′

end
@recipe function f(pq::AbstractRationalFunction{T}, a=nothing, b=nothing) where {T <: Real}
    rational_function_trim(pq, a, b)
end

