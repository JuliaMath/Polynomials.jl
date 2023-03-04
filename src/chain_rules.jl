import ChainRulesCore

function ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
    (_, ṗ, ẋ),
    ::AbstractPolynomial,
    p,
    x,

)
    @show :hi
    p(x), derivative(p)(x)
end

# ## modified from
# ## https://github.com/gdalle/ImplicitDifferentiation.jl/blob/main/src/implicit_function.jl
# function ChainRulesCore.rrule(
#     rc::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasReverseMode},
#     ::typeof(solve),
#     ZP::ZeroProblem,
#     M::AbstractUnivariateZeroMethod,
#     p;
#     kwargs...,
# )
#     xᵅ = solve(ZP, M, p; kwargs...)

#     f(x, p) = first(Callable_Function(M, ZP.F, p)(x))
#     _, pullback_f = ChainRulesCore.rrule_via_ad(rc, f, xᵅ, p)
#     _, fx, fp = pullback_f(true)
#     yp = -fp / fx

#     function pullback_solve_ZeroProblem(dy)
#         dp = yp * dy
#         return (
#             ChainRulesCore.NoTangent(),
#             ChainRulesCore.NoTangent(),
#             ChainRulesCore.NoTangent(),
#             dp,
#         )
#     end

#     return xᵅ, pullback_solve_ZeroProblem
# end
