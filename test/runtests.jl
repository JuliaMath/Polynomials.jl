# assert file to test polynomial implementation
using Compat
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
X = Poly([0.0, 1.0])
p1000 = Poly(randn(1000))

@test length(pNULL) == 1
@test length(p1000-p1000) == 1
@test length(p1000^0) == 1
@test length(0*p1000) == 1
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

@test map(degree, [pNULL,p0,p1,p2,p3,p4,p5,pN,pR,p1000]) == [0,0,0,1,2,3,4,4,2,999]

@test polyval(pN, -.125) == 276.9609375
@test polyval(pNULL, 10) == 0
@test polyval(p0, -10) == 0
@test isa(polyval(p0, X), Poly)
@test isa(polyval(p1, X), Poly)
@test polyval(poly([1//2, 3//2]), 1//2) == 0//1
@test polyder(polyint(pN)) == pN
@test polyder(pR) == Poly([-2//1,2//1])
@test polyint(pNULL,1) == p1
@test polyint(Poly(Rational[1,2,3])) == Poly(Rational[0, 1, 1, 1])
@test polyder(p3) == Poly([2,2])
@test polyder(p1) == polyder(p0) == polyder(pNULL) == pNULL

if VERSION >= v"0.4"
    @test pN(-.125) == 276.9609375
    @test pN([0.1, 0.2, 0.3]) == polyval(pN, [0.1, 0.2, 0.3])
end

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

# Test conjugation
@test conj(Poly([1im])) == Poly([-1im])
@test conj(Poly([2.0, 1im, 1+2im])) == Poly([2.0, -1im, 1-2im])
@test conj(Poly([1, 2])) == Poly([1, 2])


p0 = Poly([0])
p1 = Poly([1])
p2 = Poly([5, 6, -3, 2 ,4])
p3 = Poly([7, -3, 2, 6])
p4 = p2 * p3
@test divrem(p4, p2) == (p3, zero(p3))
@test p3%p2 == p3
@test all((abs((p2 ÷ p3 - Poly([1/9,2/3])).a)) .< eps())
@test divrem(p0,p1) == (p0,p0)
@test divrem(p1,p1) == (p1,p0)
@test divrem(p2,p2) == (p1,p0)
@test divrem(pR, pR) == (one(pR), zero(pR))
@test_throws DivideError p1 ÷ p0
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
@test_throws ErrorException pS1 ÷ pX
@test_throws ErrorException pS1 % pX

#Testing copying.
pcpy1 = Poly([1,2,3,4,5], :y)
pcpy2 = copy(pcpy1)
@test pcpy1 == pcpy2

#Tests for Pade approximants

println("Test for the exponential function.")
a = Poly(1.//convert(Vector{BigInt},gamma(BigFloat(1):BigFloat(17))),"x")
PQexp = Pade(a,8,8)
@test isapprox(convert(Float64, padeval(PQexp,1.0)), exp(1.0))
@test isapprox(convert(Float64, padeval(PQexp,-1.0)), exp(-1.0))

println("Test for the sine function.")
b = Poly(convert(Vector{BigInt},sinpi((0:16)/2)).//convert(Vector{BigInt},gamma(BigFloat(1):BigFloat(17))),"x")
PQsin = Pade(b,8,7)
@test isapprox(convert(Float64, padeval(PQsin,1.0)), sin(1.0))
@test isapprox(convert(Float64, padeval(PQsin,-1.0)),sin(-1.0))

println("Test for the cosine function.")
c = Poly(convert(Vector{BigInt},sinpi((1:17)/2)).//convert(Vector{BigInt},gamma(BigFloat(1):BigFloat(17))),"x")
PQcos = Pade(c,8,8)
@test isapprox(convert(Float64, padeval(PQcos,1.0)), cos(1.0))
@test isapprox(convert(Float64, padeval(PQcos,-1.0)), cos(-1.0))

println("Test for the summation of a factorially divergent series.")
d = Poly(convert(Vector{BigInt},(-1).^(0:60).*gamma(BigFloat(1):BigFloat(61.0))).//1,"x")
PQexpint = Pade(d,30,30)
@compat println("The approximate sum of the divergent series is:  ", Float64(padeval(PQexpint,1.0)))
println("The approximate sum of the convergent series is: ",exp(1)*(-γ-sum([(-1).^k/k./gamma(k+1) for k=1:20])))
@test isapprox(convert(Float64, padeval(PQexpint,1.0)),
               exp(1)*(-γ-sum([(-1).^k/k./gamma(k+1) for k=1:20])))


## polyfit
xs = linspace(0, pi, 10)
ys = sin(xs)
p = polyfit(xs, ys)
p = polyfit(xs, ys, :t)
p = polyfit(xs, ys, 2)
@test maximum(abs(map(x->polyval(p, x), xs) - ys)) <= 0.03


## truncation
p1 = Poly([1,1]/10)
p2 = Poly([1,2]/10)
p3 = Poly([1,3]/10)
psum = p1 + p2 - p3
@test degree(psum) == 1         # will have wrong degree
@test degree(truncate(psum)) == 0 # the degree should be correct after truncation

@test truncate(Poly([2,1]),reltol=1/2,abstol=0) == Poly([2])
@test truncate(Poly([2,1]),reltol=1,abstol=0)   == Poly([0])
@test truncate(Poly([2,1]),reltol=0,abstol=1)   == Poly([2])


## setindex!
println("Test for setindex!()")
p1    = Poly([1,2,1])
p1[5] = 1
@test p1[5] == 1
@test p1 == Poly([1,2,1,0,0,1])


## elementwise operations #52
println("Test for element-wise operations")
p1  = Poly([1, 2])
p2  = Poly([3, 1.])
p   = [p1, p2]
q   = [3, p1]
@test isa(q,Vector{Poly{Int64}})
psum  = p+3
pprod = p*3
pmin  = p-3
@test isa(psum, Vector{Poly{Float64}})
@test isa(pprod,Vector{Poly{Float64}})
@test isa(pmin, Vector{Poly{Float64}})

## getindex with ranges #43
p1 = Poly([4,5,6])
@test all(p1[0:1] .== [4,5])
@test all(p1[0:end] .== [4,5,6])
p1[0:1] = [7,8]
@test all(p1[0:end] .== [7,8,6])
