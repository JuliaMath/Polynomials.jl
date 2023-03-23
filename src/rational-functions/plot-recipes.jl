## Plot recipe
## define a hueristic to work around asymptotes
## just sort of successful
@recipe function f(pq::AbstractRationalFunction{T}, a=nothing, b=nothing) where {T}

    xlims = get(plotattributes, :xlims, (nothing, nothing))
    ylims = get(plotattributes, :ylims, (nothing, nothing))
    rational_function_trim(pq, a, b, xlims, ylims)    

end

isapproxreal(x::Real) = true
isapproxreal(x::Complex{T}) where {T} = imag(x) <= sqrt(eps(real(T)))
function toobig(pq)
    x -> begin
        y = pq(x)
        isinf(y) && return true
        isnan(y) && return true
        abs(y) > 1e8 && return true
        return false
    end
end

function rational_function_trim(pq, a, b, xlims, ylims)

    p,q = lowest_terms(//(pq...), method=:numerical)
    dpq = derivative(p//q)
    p′,q′ = lowest_terms(dpq)

    λs = Multroot.multroot(q).values
    λs = isempty(λs) ? λs : real.(filter(isapproxreal, λs))

    cps = Multroot.multroot(p′).values
    cps = isempty(cps) ? cps : real.(filter(isapproxreal, cps))
    cps = isempty(cps) ? cps : filter(!toobig(pq), cps)

    a = isnothing(a) ? xlims[1] : a
    b = isnothing(b) ? xlims[2] : b
    
    if isnothing(a) && isnothing(b)
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

    n = 601
    xs = range(a,stop=b, length=n)
    ys = pq.(xs)
    Mcps = isempty(cps) ? 5 : 3*maximum(abs, pq.(cps))
    M = max(5, Mcps, 1.25*maximum(abs, pq.((a,b))))

    lo = isnothing(ylims[1]) ? -M : ylims[1]
    hi = isnothing(ylims[2]) ?  M : ylims[2]
    ys′ = [lo <= yᵢ <= hi ? yᵢ : NaN for yᵢ ∈ ys]
    xs, ys′

end

