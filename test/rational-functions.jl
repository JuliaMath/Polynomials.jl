using Polynomials
using Test
using LinearAlgebra

@testset "Basic constructor, arithmetic" begin

    p,q = fromroots(Polynomial, [1,2,3]), fromroots(Polynomial, [2,3,4])
    r,s = SparsePolynomial([1,2,3], :x), SparsePolynomial([1,2,3], :y)
    t,u = Polynomial([1,2,pi]), SparsePolynomial([1.1, 2.2, 3.4], :y)
    
    # constructor
    @test p // q isa RationalFunction
    @test p // r isa RationalFunction
    @test_throws ArgumentError r // s
    @test RationalFunction(p) == p // one(p)

    # We expect p::P // q::P (same type polynomial).
    # As Immutable Polynomials have N as type parameter, we disallow
    @test_throws ArgumentError variable(ImmutablePolynomial) // variable(ImmutablePolynomial)
    
    pq = p // t # promotes to type of t
    @test pq isa RationalFunction{Float64, :x}
   
    # iterations, broadcast
    pp,qq = p // q
    @test (pp,qq) == (p,q)
    @test eltype(p//q) == typeof(pp)
    u = gcd(p//q...)
    @test u/u[end] ≈ fromroots(Polynomial, [2,3])
    @test degree.(p//q) == degree(p//q) # no broadcast over rational functions
    
    # evaluation
    pq = p//q
    @test pq(2.5) ≈ p(2.5) / q(2.5)
    @test pq(2) ≈ fromroots([1,3])(2) / fromroots([3,4])(2)
    
    # arithmetic
    rs = r // (r-1)
    x = 1.2
    @test (pq + rs)(x) ≈ (pq(x) + rs(x))
    @test (pq * rs)(x) ≈ (pq(x) * rs(x))
    @test (-pq)(x) ≈ -p(x)/q(x)
    @test pq .* (pq, pq) == (pq*pq, pq*pq)
    @test pq .* [pq pq] == [pq*pq pq*pq]
    
    
    # derivative
    pq = p // one(p)
    x = 1.23
    @test derivative(pq)(x) ≈ derivative(p)(x)
    pq = p // q
    dp,dq = derivative.((p,q))
    @test derivative(pq)(x) ≈ (dp(x)*q(x) - p(x)*dq(x))/q(x)^2
    
    # integral
    pq = p // (2*one(p))
    @test iszero(derivative(integrate(pq)) - pq)
    pq = p // q
    @test_throws ArgumentError integrate(pq) # need logs terms in general
    
    # lowest terms
    pq = p // q
    pp, qq = Polynomials.pqs(lowest_terms(pq))
    @test all(abs.(pp.(roots(qq))) .> 1/2)

    # ≈
    @test (2p)//(2one(p)) ≈ p//one(p)
    @test p//one(p) ≈  p//one(p) + eps()
    x = variable(p)
    @test (x*p)//x ≈ p // one(p)

    
end

@testset "zeros, poles, residues" begin
    x = variable(Polynomial{Float64})
    p,q = x^7,(x^5 + 4*x^4 + 7*x^3 + 8*x^2 + 5*x + 2)
    pq = p // q


    ## zeros
    zs = roots(pq)
    @test length(zs.zs) == 1
    @test zs.multiplicities == [7]
    
    ## poles
    ps = poles(pq)
    @test length(ps.zs) == 3
    
    ## residues
    d,r = residues(pq)
    @test d ≈ 9 - 4x + x^2

    z = variable(pq)
    for (λ, rs) ∈ r # reconstruct p/q from output of `residues`
        for (i,rᵢ) ∈ enumerate(rs)
            d += rᵢ/(z-λ)^i
        end
    end
    @test norm(numerator(lowest_terms(d - pq)), Inf) <= sqrt(eps())

end

@testset "As matrix elements" begin
    
    p, q = Polynomial([1,2], :x), Polynomial([1,2],:y)
    pp = p // (p-1)
    PP = typeof(pp)
    
    r, s = SparsePolynomial([1,2], :x), SparsePolynomial([1,2],:y)
    rr = r // (r-1)

    ## T, Polynomial{T} promotes
    @test eltype([1, p, pp]) == PP

    ## test mixed types promote polynomial type
    @test eltype([pp rr p r]) == PP

    ## test non-constant polys of different symbols throw error
    @test_throws ArgumentError [pp, q]
    @test_throws ArgumentError [pp, q//(q-1)]

    ## constant polys will work -- if specified
    @test eltype(PP[pp one(q)]) == PP
    @test eltype(PP[pp one(q//(q-1))]) == PP

    ## if not specified, the rational function will not promote
    @test eltype([pp one(q)]) == PP
    @test_throws ArgumentError [pp one(q//(q-1))]

end


@testset "Rational function fit" begin
    
    p,q = Polynomial([1,1,1]), Polynomial([1,2])
    xs = range(-0.25,stop=1, length=15)
    ys = (p//q).(xs)
    pq = fit(RationalFunction, xs, ys, 2, 2)
    @test numerator(pq) ≈ p
    @test denominator(pq) ≈ q

    x = variable(Polynomial{Float64})
    pq = (1+x)//(1-x)
    xs = 2.0:.1:3;
    ys = pq.(xs);
    v = fit(RationalFunction, xs, ys, 2, 2)
    @test maximum(abs, v(x)-pq(x) for x ∈ 2.1:0.1:3.0) <= sqrt(eps())
    
end

@testset "Pade fit" begin

    #Tests for Pade approximants from Pade.jl

    # exponential
    a = Polynomial( 1 .// factorial.(big(0):17) )
    pq = fit(RationalFunction,a,8,8)
    @test pq(1.0) ≈ exp(1.0)
    @test pq(-1.0) ≈ exp(-1.0)    

    # sine
    b = Polynomial(Int.(sinpi.((0:16)/2)) .// factorial.(big(0):16))
    pq = fit(RationalFunction, b, 7, 8)
    @test pq(1.0) ≈ sin(1.0)
    @test pq(-1.0) ≈ sin(-1.0)    

    # cosine
    c = Polynomial(Int.(sinpi.((1:17)/2)) .// factorial.(big(0):16))
    pq = fit(RationalFunction, c, 8, 8)
    @test pq(1.0) ≈ cos(1.0)
    @test pq(-1.0) ≈ cos(-1.0)    

    # summation of a factorially divergent series
    d = Polynomial( (-big(1)).^(0:60) .* factorial.(big(0):60) )
    pq = fit(RationalFunction, d, 30, 30)
    @test pq(1.0) ≈ exp(1)*(-Base.MathConstants.γ-sum([(-1).^k/k./factorial(k) for k=1:20]))

end
