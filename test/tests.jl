# assert file to test polynomial implementation
using Base.Test
using Polynomials

pNULL = Poly(Float32[])
p0 = Poly([0])
p1 = Poly([1,0,0,0,0,0,0,0,0,0,0,0,0,0])
p2 = Poly([1,1,0,0])
p3 = Poly([1,2,1,0,0,0,0])
p4 = Poly([1,3,3,1,0,0])
p5 = Poly([1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
pN = Poly([276,3,87,15,24,0])
pR = Poly([3//4, -2//1, 1//1])
p1000 = Poly(randn(1000))

@test length(pNULL) == 0
@test length(p1000) == 1000
sprint(show, p1000)
sprint(show, pNULL)

@test p3 == Poly([1,2,1])
@test pN*10 == Poly([2760, 30, 870, 150, 240])
@test pN/10 == Poly([27.6, 0.3, 8.7, 1.5, 2.4])
@test 10*pNULL + pN == pN
@test 10*p0 + pN == pN
@test p5 + 2*p1 == Poly([3,4,6,4,1])
@test 10*pNULL - pN == -pN
@test p0 - pN == -pN
@test p5 - 2*p1 == Poly([-1,4,6,4,1])
@test p2*p2*p2 == p4
@test p2^4 == p5
@test pNULL^3 == pNULL
@test pNULL*pNULL == pNULL

@test map(Polynomials.degree, [pNULL,p0,p1,p2,p3,p4,p5,pN,pR,p1000]) == [0,0,0,1,2,3,4,4,2,999] 

@test polyval(pN, -.125) == 276.9609375
@test polyval(pNULL, 10) == 0
@test polyval(p0, -10) == 0
@test polyval(poly([1//2, 3//2]), 1//2) == 0//1
@test polyder(polyint(pN)) == pN
@test polyder(pR) == Poly([-2//1,2//1])
@test polyint(pNULL,1) == p1
@test polyint(Poly(Rational[1,2,3])) == Poly(Rational[0, 1, 1, 1])
@test polyder(p3) == Poly([2,2])
@test polyder(p1) == polyder(p0) == polyder(pNULL) == pNULL

@test poly([-1,-1]) == p3
@test roots(p0)==roots(p1)==roots(pNULL)==[] 
@test roots(p2) == [-1]
a_roots = copy(pN.a)
@test all(abs(sort(roots(poly(a_roots))) - sort(a_roots)) .< 1e6)
@test length(roots(p5)) == 4
@test roots(pNULL) == []
@test roots(pR) == [1//2, 3//2]

@test pNULL + 2 == p0 + 2 == 2 + p0 == Poly([2])
@test p2 - 2 == -2 + p2 == Poly([-1,1])
@test 2 - p2 == Poly([1,-1])

p0 = Poly([0])
p1 = Poly([1])
p2 = Poly([5, 6, -3, 2 ,4])
p3 = Poly([7, -3, 2, 6])
p4 = p2 * p3
@test divrem(p4, p2) == (p3, zero(p3))
@test p3%p2 == p3
@test all((abs((p2/p3 - Poly([1/9,2/3])).a)) .< eps())
@test divrem(p0,p1) == (p0,p0)
@test divrem(p1,p1) == (p1,p0)
@test divrem(p2,p2) == (p1,p0)
@test divrem(pR, pR) == (one(pR), zero(pR))
@test_throws DivideError p1/p0
@test_throws DivideError divrem(p0,p0)

#Tests for multivariable support
pX = Poly([1, 2, 3, 4, 5])
pS1 = Poly([1, 2, 3, 4, 5], "s")
pS2 = Poly([1, 2, 3, 4, 5], 's')
pS3 = Poly([1, 2, 3, 4, 5], :s)
@test pX != pS1
@test pS1 == pS2
@test pS1 == pS3
@test_throws ErrorException pS1 + pX
@test_throws ErrorException pS1 - pX
@test_throws ErrorException pS1 * pX
@test_throws ErrorException pS1 / pX
@test_throws ErrorException pS1 % pX

#Tests for Pade approximants

println("Test for the exponential function.")
a = Poly(1.//convert(Vector{Int},gamma(1:17)),"x")
PQexp = Pade(a,8,8)
@test padeval(PQexp,1.0) == exp(1.0)
@test padeval(PQexp,-1.0) == exp(-1.0)

println("Test for the sine function.")
b = Poly(convert(Vector{BigInt},sinpi((0:16)/2)).//convert(Vector{BigInt},gamma(BigFloat("1.0"):BigFloat("17.0"))),"x")
PQsin = Pade(b,8,7)
@test isapprox(padeval(PQsin,1.0),sin(1.0))
@test isapprox(padeval(PQsin,-1.0),sin(-1.0))

println("Test for the cosine function.")
c = Poly(convert(Vector{BigInt},sinpi((1:17)/2)).//convert(Vector{BigInt},gamma(BigFloat("1.0"):BigFloat("17.0"))),"x")
PQcos = Pade(c,8,8)
@test isapprox(padeval(PQcos,1.0),cos(1.0))
@test isapprox(padeval(PQcos,-1.0),cos(-1.0))

println("Test for the summation of a factorially divergent series.")
d = Poly(convert(Vector{BigInt},(-1).^(0:60).*gamma(BigFloat("1.0"):BigFloat("61.0"))).//1,"x")
PQexpint = Pade(d,30,30)
println("The approximate sum of the divergent series is:  ",float64(padeval(PQexpint,1.0)))
println("The approximate sum of the convergent series is: ",exp(1)*(-γ-sum([(-1).^k/k./gamma(k+1) for k=1:20])))
@test isapprox(padeval(PQexpint,1.0) , exp(1)*(-γ-sum([(-1).^k/k./gamma(k+1) for k=1:20])))

#Tests for orthogonal polynomials

function orthogonal_poly_test(F,Ps,num=[],test_small=false)
    if num == []
        num = ones(length(Ps))
    end
    for n=0:length(Ps)-1
        diff = F(n) - Poly(Ps[n+1])/num[n+1]
        # println("$(n): $(F(n)) - ($(Poly(Ps[n+1])/num[n+1])) = $(diff) @ 1.1 = $(polyval(diff, 1.1))")
        if test_small == false
            @test F(n) == Poly(Ps[n+1])/num[n+1]
        else
            # Since some fractions cannot be represented exactly by
            # floating point arithmetic, we test if the polynomials are
            # "close" at a certain point of evaluation.
            @test_approx_eq_eps polyval(diff, test_small) 0 1e-14
        end
    end
end

println("Testing the first 11 Legendre polynomials")
orthogonal_poly_test(Orthogonal.P, {[1] [0,1] [-1,0,3] [0,-3,0,5] [3,0,-30,0,35] [0,15,0,-70,0,63] [-5,0,105,0,-315,0,231] [0,-35,0,315,0,-693,0,429] [35,0,-1260,0,6930,0,-12012,0,6435] [0,315,0,-4620,0,18018,0,-25740,0,12155] [-63,0,3465,0,-30030,0,90090,0,-109395,0,46189]}, [1 1 2 2 8 8 16 16 128 128 256])

println("Testing the first 16 Chebyshev polynomials of the first kind")
orthogonal_poly_test(Orthogonal.T, {[1] [0,1] [-1,0,2] [0,-3,0,4] [1,0,-8,0,8] [0,5,0,-20,0,16] [-1,0,18,0,-48,0,32] [0,-7,0,56,0,-112,0,64] [1,0,-32,0,160,0,-256,0,128] [0,9,0,-120,0,432,0,-576,0,256] [-1,0,50,0,-400,0,1120,0,-1280,0,512] [0,-11,0,220,0,-1232,0,2816,0,-2816,0,1024] [1,0,-72,0,840,0,-3584,0,6912,0,-6144,0,2048] [0,13,0,-364,0,2912,0,-9984,0,16640,0,-13312,0,4096] [-1,0,98,0,-1568,0,9408,0,-26880,0,39424,0,-28672,0,8192] [0,-15,0,560,0,-6048,0,28800,0,-70400,0,92160,0,-61440,0,16384]})

println("Testing the first 16 Chebyshev polynomials of the second kind")
orthogonal_poly_test(Orthogonal.U, {[1] [0,2] [-1,0,4] [0,-4,0,8] [1,0,-12,0,16] [0,6,0,-32,0,32] [-1,0,24,0,-80,0,64] [0,-8,0,80,0,-192,0,128] [1,0,-40,0,240,0,-448,0,256] [0,10,0,-160,0,672,0,-1024,0,512] [-1,0,60,0,-560,0,1792,0,-2304,0,1024] [0,-12,0,280,0,-1792,0,4608,0,-5120,0,2048] [1,0,-84,0,1120,0,-5376,0,11520,0,-11264,0,4096]})

println("Testing the first 11 Hermite polynomials")
orthogonal_poly_test(Orthogonal.He, {[1] [0,1] [-1,0,1] [0,-3,0,1] [3,0,-6,0,1] [0,15,0,-10,0,1] [-15,0,45,0,-15,0,1] [0,-105,0,105,0,-21,0,1] [105,0,-420,0,210,0,-28,0,1] [0,945,0,-1260,0,378,0,-36,0,1] [-945,0,4725,0,-3150,0,630,0,-45,0,1]})
orthogonal_poly_test(Orthogonal.H, {[1] [0,2] [-2,0,4] [0,-12,0,8] [12,0,-48,0,16] [0,120,0,-160,0,32] [-120,0,720,0,-480,0,64] [0,-1680,0,3360,0,-1344,0,128] [1680,0,-13440,0,13440,0,-3584,0,256] [0,30240,0,-80640,0,48384,0,-9216,0,512] [-30240,0,302400,0,-403200,0,161280,0,-23040,0,1024]})

println("Testing the first 11 Laguerre polynomials")
orthogonal_poly_test(Orthogonal.L, {[1] [1,-1] [2,-4,1] [6,-18,9,-1] [24,-96,72,-16,1] [120,-600,600,-200,25,-1] [720,-4320,5400,-2400,450,-36,1] [5040,-35280,52920,-29400,7350,-882,49,-1] [40320,-322560,564480,-376320,117600,-18816,1568,-64,1]}, [1 1 2 6 24 120 720 5040 40320], 1.1)
