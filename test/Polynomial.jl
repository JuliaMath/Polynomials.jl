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
@test polydir(polyint(pN)) == pN
@test polyint(pNULL,1) == p1
@test polydir(p3) == Poly([2,2])
@test polydir(p1) == polydir(p0) == polydir(pNULL) == pNULL

@test poly([-1,-1]) == p3
@test roots(p0)==roots(p1)==roots(pNULL)==[] 
@test roots(p2) == [-1]
a_roots = copy(pN.a)
@test all(abs(sort(roots(poly(a_roots))) - sort(a_roots)) .< 1e6)
@test length(roots(p5)) == 4
@test roots(pNULL) == []

@test pNULL + 2 == p0 + 2 == 2 + p0 == Poly([2])
@test p2 - 2 == -2 + p2 == Poly([1,-1])
@test 2 - p2 == Poly([-1,1])

