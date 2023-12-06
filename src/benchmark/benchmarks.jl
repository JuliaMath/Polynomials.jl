using BenchmarkTools
using Polynomials

Ps = (P=Polynomial, IP=ImmutablePolynomial, SP=SparsePolynomial,
      LP=LaurentPolynomial)

const SUITE = BenchmarkGroup()

for P ∈ Ps
    P′ = string(P)
    SUITE[P′] = BenchmarkGroup(["Polynomials", P′])
    SUITE[P′]["evaluation"] =
        (@benchmarkable p(c) setup=(p=$P(rand(20)); c=rand()))
    SUITE[P′]["scalar addition"] =
        (@benchmarkable p + c setup=(p=$P(rand(20)); c=rand()))
    SUITE[P′]["scalar multiplication"] =
        (@benchmarkable p * c setup=(p=$P(rand(20)); c=rand()))
    SUITE[P′]["scalar divistion"] =
        (@benchmarkable p / c setup=(p=$P(rand(20)); c=rand()))
    SUITE[P′]["poly addition"] =
        (@benchmarkable p + q setup=(p=$P(rand(20)); q = $P(rand(25))))
    SUITE[P′]["poly multiplication"] =
        (@benchmarkable p * q setup=(p=$P(rand(20)); q = $P(rand(25))))
    SUITE[P′]["mixed var add"] =
        (@benchmarkable p + q setup=(p=$P(rand(20)); q=$P(rand(1),:y)))
    SUITE[P′]["mixed var mul"] =
        (@benchmarkable p * q setup=(p=$P(rand(20)); q=$P(rand(1),:y)))
    SUITE[P′]["differentiation"] =
        (@benchmarkable derivative(p) setup=(p=$P(rand(20))))
    SUITE[P′]["integration"] =
        (@benchmarkable integrate(p) setup=(p=$P(rand(20))))
end
