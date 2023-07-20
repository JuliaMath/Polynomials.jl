## XXX move to contrib.jl
gtτ(x, τ) = abs(x) > τ

# return index or nothing of last non "zdero"
# chop! then considers cases i==nothing, i=length(x), i < length(x)
function chop_right_index(x; rtol=nothing, atol=nothing)

    isempty(x) && return nothing
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, norm(x,2) * δ)
    i = findlast(Base.Fix2(gtτ, τ), x)
    i

end

function chop_left_index(x; rtol=nothing, atol=nothing)
    isempty(x) && return nothing
    δ = something(rtol,0)
    ϵ = something(atol,0)
    τ = max(ϵ, norm(x,2) * δ)
    i = findfirst(Base.Fix2(gtτ,τ), x)
    i
end

XXX() = throw(ArgumentError("Method not defined"))
