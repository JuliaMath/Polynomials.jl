module PolynomialsChainRulesCoreExt

using Polynomials
import ChainRulesCore

function ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
    (_, Δx),
    p::AbstractPolynomial,
    x
)
    p(x), derivative(p)(x)*Δx
end

function ChainRulesCore.rrule(p::AbstractPolynomial, x)
    _pullback(ΔΩ) = (ChainRulesCore.NoTangent(), derivative(p)(x))
    return (p(x), _pullback)
end

end
