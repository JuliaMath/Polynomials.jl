module Orthogonal
using Polynomials

macro recur_poly_gen(name, rule, p₁)
    generator = gensym()
    @eval begin
        ($name)(n::Integer, x = Poly([0, 1])) = ($generator)(n::Integer, x)[1]
        function ($generator){T}(n::Integer, x::T)
            if n == 0
                one(x), zero(x)
            elseif n == 1
                $p₁, one(x)
            else
                p₋₁, p₋₂ = ($generator)(n-1, x)
                $rule, p₋₁
            end
        end
    end
end

# Jacobi polynomials are orthogonal on the interval [-1,1] with respect
# to the weight (1-x)^α(1+x)^β. Legendre and Chebyshev polynomials are
# special cases.
# TODO: General Jacobi, Gegenbauer and associated Legendre polynomials.
@recur_poly_gen(T, 2x*p₋₁ - p₋₂, x) # Chebyshev polynomials, first kind
@recur_poly_gen(U, 2x*p₋₁ - p₋₂, 2x) # Chebyshev polynomials, second kind
@recur_poly_gen(P, ((2n - 1)x*p₋₁ - (n - 1)p₋₂)//n, x) # Legendre polynomials

# Hermite polynomials are orthogonal on the interval (-∞,∞) with
# weight exp(-x²/a) where a = 2 for probabilists and a = 1 for
# physicists.
@recur_poly_gen(He, x*p₋₁ - (n-1)p₋₂, x) # Probabilists' Hermite polynomials
@recur_poly_gen(H, 2x*p₋₁ - 2(n-1)p₋₂, 2x) # Physicists' Hermite polynomials

# Associated Laguerre polynomials are orthogonal on the interval [0,∞)
# with respect to the weight x^α exp(-x).
# TODO: Associated Laguerre polynomials.
@recur_poly_gen(L, ((2n - 1 - x)p₋₁ - (n-1)p₋₂)//n, 1-x)
end
