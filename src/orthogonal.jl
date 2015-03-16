module Orthogonal
using Polynomials
#=
Jacobi polynomials are orthogonal on the interval [-1,1] with respect
to the weight (1-x)^α(1+x)^β. Legendre and Chebyshev polynomials are
special cases.

They are defined through the recursion relation:

2n(n + α + β)(2n + α + β - 2)Pₙ{α,β}(x) =
    (2n + α + β - 1)[(2n + α + β)(2n + α + β - 2)x + α² - β²]Pₙ₋₁{α,β}(x)
    - 2(n + α - 1)(n + β - 1)(2n + α + β)Pₙ₋₂{α,β}(x)

with

P₀{α,β}(x) = 1

and

P₁{α,β} = 1/2[2(α + 1) + (α + β + 2)(x - 1)]


 =#

function P(n::Integer, α, β)
    n == 0 ? Poly([1]) : n == 1 ? (α + 1) + 1/2*(α + β + 2)*Poly([-1,1]) :
    ((2n + α + β - 1)*((2n + α + β)*(2n + α + β - 2)*Poly([0,1]) + α^2 - β^2)*P(n-1,α,β)
     - 2(n + α - 1)*(n + β - 1)*(2n + α + β)*P(n-2,α,β))/(2n*(n + α + β)*(2n + α + β - 2))
end

# Legendre polynomials
function P(n::Integer)
    P(n,0,0)
end

# Chebyshev polynomials of the first kind
T = n::Integer -> n == 0 ? Poly([1]) : n == 1 ? Poly([0,1]) : Poly([0,2])*T(n-1) - T(n-2)

# Chebyshev polynomials of the second kind
U = n::Integer -> n == 0 ? Poly([1]) : n == 1 ? Poly([0,2]) : Poly([0,2])*U(n-1) - U(n-2)

####################################################################################

# Hermite polynomials are orthogonal on the interval (-∞,∞) with
# weight exp(-x²/a) where a = 2 for probabilists and a = 1 for
# physicists.

# Probabilists' Hermite polynomials:
# Heₙ₊₁(x) = xHeₙ(x) - nHeₙ₋₁(x)
function He(n::Integer)
    n == 0 ? Poly([1]) : n == 1 ? Poly([0, 1]) : Poly([0,1])*He(n-1) - (n-1)*He(n-2)
end

# Physicists' Hermite polynomials:
# Hₙ₊₁(x) = 2xHₙ(x) - 2nHₙ₋₁(x)
function H(n::Integer)
    n == 0 ? Poly([1]) : n == 1 ? Poly([0, 2]) : Poly([0, 2])*H(n-1) - 2(n-1)*H(n-2)
end

####################################################################################

#=
Associated Laguerre polynomials are orthogonal on the interval [0,∞)
with respect to the weight x^α exp(-x).

They are defined through the recursion relation:

Lₙ₊₁{α}(x) = [(2n + 1 + α - x)Lₙ{α}(x) - (n + α)Lₙ₋₁{α}(x)]/(n + 1)

with

L₀{α}(x) = 1

and

L₁{α}(x) = 1 + α - x
=#

function L(n::Integer, α)
    n == 0 ? Poly([1]) : n == 1 ? Poly([1 + α, -1]) : (Poly([2n - 1 + α, -1])*L(n-1, α) - (n - 1 + α)*L(n-2, α))/n
end

# Laguerre polynomials
function L(n::Integer)
    L(n, 0)
end
end
