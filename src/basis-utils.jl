_norm(x,p=2) = real(sqrt(sum(xᵢ^2 for xᵢ ∈ x)))
gtτ(x, τ) = abs(x) > τ


# return index or nothing of last non "zdero"
# chop! then considers cases i==nothing, i=length(x), i < length(x)
function chop_right_index(x; rtol=nothing, atol=nothing)

    isempty(x) && return nothing
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, _norm(x,2) * δ)
    i = findlast(Base.Fix2(gtτ, τ), x)
    i

end



function chop_left_index(x; rtol=nothing, atol=nothing)
    isempty(x) && return nothing
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, _norm(x,2) * δ)
    i = findfirst(Base.Fix2(gtτ,τ), x)
    i
end

## put here, not with type definition, in case reuse is possible
## `conv` can be used with matrix entries, unlike `fastconv`
function XXXconv(p::Vector{T}, q::Vector{S}) where {T,S}
    (isempty(p) || isempty(q)) && return promote_type(T, S)[]
    as = [p[1]*q[1]]
    z = zero(as[1])
    n,m = length(p)-1, length(q)-1
    for i ∈ 1:n+m
        Σ = z
        for j ∈ max(0, i-m):min(i,n)
            Σ = muladd(p[1+j], q[1 + i-j], Σ)
        end
        push!(as, Σ)
    end
    as
end

XXX() = throw(ArgumentError("Method not defined"))
