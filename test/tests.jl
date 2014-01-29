# assert file to test polynomial implementation
using Base.Test
using Polynomial

pNULL = Poly(Float32[])
p0 = Poly([0])
p1 = Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,1])
p2 = Poly([0,0,1,1])
p3 = Poly([0,0,0,0,1,2,1])
p4 = Poly([0,1,3,3,1])
p5 = Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,6,4,1])
pN = Poly([0,24,15,87,3,276])
pR = Poly([1//1, -2//1, 3//4])
p1000 = Poly(randn(1000))

@test length(pNULL) == 0
@test length(p1000) == 1000
sprint(show, p1000)
sprint(show, pNULL)

@test p3 == Poly([1,2,1])
@test pN*10 == Poly([240, 150, 870, 30, 2760])
@test pN/10 == Poly([2.4, 1.5, 8.7, 0.3, 27.6])
@test 10*pNULL + pN == pN
@test 10*p0 + pN == pN
@test p5 + 2*p1 == Poly([1,4,6,4,3])
@test 10*pNULL - pN == -pN
@test p0 - pN == -pN
@test p5 - 2*p1 == Poly([1,4,6,4,-1])
@test p2*p2*p2 == p4
@test p2^4 == p5
@test pNULL^3 == pNULL
@test pNULL*pNULL == pNULL

@test polyval(pN, -.125) == 276.9609375
@test polyval(pNULL, 10) == 0
@test polyval(p0, -10) == 0
@test polyval(poly([1//2, 3//2]), 1//2) == 0//1
@test polyder(polyint(pN)) == pN
@test polyder(pR) == Poly([2//1, -2//1])
@test polyint(pNULL,1) == p1
@test polyint(Poly(Rational[3, 2, 1])) == Poly(Rational[1, 1, 1, 0])
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
@test p2 - 2 == -2 + p2 == Poly([1,-1])
@test 2 - p2 == Poly([-1,1])

p0 = Poly([0])
p1 = Poly([1])
p2 = Poly([4, 2, -3, 6, 5])
p3 = Poly([6, 2, -3, 7])
p4 = p2 * p3
@test divrem(p4, p2) == (p3, zero(p3))
@test p3%p2 == p3
@test all((abs((p2/p3 - Poly([2/3,1/9])).a)) .< eps())
@test divrem(p0,p1) == (p0,p0)
@test divrem(p1,p1) == (p1,p0)
@test divrem(p2,p2) == (p1,p0)
@test divrem(pR, pR) == (one(pR), zero(pR))
@test_throws p1/p0
@test_throws divrem(p0,p0)
