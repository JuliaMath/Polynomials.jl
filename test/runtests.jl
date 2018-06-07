# assert file to test polynomial implementation
using Compat
using Compat.Test
using Compat.LinearAlgebra
using Polynomials

import Compat.SparseArrays: sparse, speye, nnz

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
T = Int64
Poly{T}([zero(T), one(T)])
Poly{T}([zero(T), one(T)], :y)

p1000 = Poly(randn(1000))

@test length(pNULL) == 1
@test length(p1000-p1000) == 1
@test length(p1000^0) == 1
@test length(0*p1000) == 1
@test length(p1000) == 1000
sprint(show, p1000)
sprint(show, pNULL)

@test Poly([1, -2, 1]) == poly(diagm(0=>[1., 1.]))

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

@test polyval(poly(Int[]), 2.) == 1.
@test polyval(pN, -.125) == 276.9609375
@test polyval(pNULL, 10) == 0
@test polyval(p0, -10) == 0
@test isa(polyval(p0, X), Poly)
@test isa(polyval(p1, X), Poly)
@test polyval(poly([1//2, 3//2]), 1//2) == 0//1
@test polyder(polyint(pN)) == pN
@test polyder(pR) == Poly([-2//1,2//1])
@test polyder(p3) == Poly([2,2])
@test polyder(p1) == polyder(p0) == polyder(pNULL) == pNULL
@test_throws ErrorException polyder(pR, -1)
@test polyint(pNULL,1) == p1
@test polyint(Poly(Rational[1,2,3])) == Poly(Rational[0, 1, 1, 1])
@test polyint(p2, 0, 2) == 4.0


@test pN(-.125) == 276.9609375
@test pN([0.1, 0.2, 0.3]) == polyval(pN, [0.1, 0.2, 0.3])

@test poly([-1,-1]) == p3
@test roots(p0)==roots(p1)==roots(pNULL)==[]
@test roots(Poly([0,1,0])) == [0.]
@test roots(p2) == [-1]
a_roots = copy(pN.a)
@test all(map(abs,sort(roots(poly(a_roots))) - sort(a_roots)) .< 1e6)
@test length(roots(p5)) == 4
@test roots(pNULL) == []
@test sort(roots(pR)) == [1//2, 3//2]
# default type in variable() is Float64
x = variable()
plarge = 8.362779449448982e41 - 2.510840694154672e57x + 4.2817430781178795e44x^2 - 1.6225927682921337e31x^3 + 1.0x^4  # #120
@test length(roots(plarge)) == 4


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
@test all((map(abs,(p2 ÷ p3 - Poly([1/9,2/3])).a)) .< eps())
@test divrem(p0,p1) == (p0,p0)
@test divrem(p1,p1) == (p1,p0)
@test divrem(p2,p2) == (p1,p0)
@test divrem(pR, pR) == (one(pR), zero(pR))
@test_throws DivideError p1 ÷ p0
@test_throws DivideError divrem(p0,p0)

# test chop
pchop = Poly([1, 2, 3, 0, 0, 0])
pchopped = chop(pchop)
@test roots(pchop) == roots(pchopped)

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
a = Poly(1 .//convert(Vector{BigInt},map(gamma,BigFloat(1):BigFloat(17))),"x")
PQexp = Pade(a,8,8)
@test isapprox(convert(Float64, padeval(PQexp,1.0)), exp(1.0))
@test isapprox(convert(Float64, padeval(PQexp,-1.0)), exp(-1.0))

println("Test for the sine function.")
b = Poly(convert(Vector{BigInt},map(sinpi,(0:16)/2)).//convert(Vector{BigInt},map(gamma,BigFloat(1):BigFloat(17))),"x")
PQsin = Pade(b,8,7)
@test isapprox(convert(Float64, padeval(PQsin,1.0)), sin(1.0))
@test isapprox(convert(Float64, padeval(PQsin,-1.0)),sin(-1.0))

println("Test for the cosine function.")
c = Poly(convert(Vector{BigInt},map(sinpi,(1:17)/2)).//convert(Vector{BigInt},map(gamma,BigFloat(1):BigFloat(17))),"x")
PQcos = Pade(c,8,8)
@test isapprox(convert(Float64, padeval(PQcos,1.0)), cos(1.0))
@test isapprox(convert(Float64, padeval(PQcos,-1.0)), cos(-1.0))

const _γ = 0.5772156649015
println("Test for the summation of a factorially divergent series.")
d = Poly(convert(Vector{BigInt},(-1).^(0:60).* map(gamma,BigFloat(1):BigFloat(61.0))).//1,"x")
PQexpint = Pade(d,30,30)
@compat println("The approximate sum of the divergent series is:  ", Float64(padeval(PQexpint,1.0)))
println("The approximate sum of the convergent series is: ",exp(1)*(-_γ-sum([(-1).^k/k./gamma(k+1) for k=1:20])))
@test isapprox(convert(Float64, padeval(PQexpint,1.0)),
               exp(1)*(-_γ-sum([(-1).^k/k./gamma(k+1) for k=1:20])))


## polyfit
if VERSION < v"0.7-"
    xs = linspace(0, pi, 10)
else
    xs = range(0, stop=pi, length=10)
end

ys = map(sin,xs)
p = polyfit(xs, ys)
p = polyfit(xs, ys, :t)
p = polyfit(xs, ys, 2)
@test maximum(map(abs,map(x->polyval(p, x), xs) - ys)) <= 0.03


## truncation
p1 = Poly([1,1]/10)
p2 = Poly([1,2]/10)
p3 = Poly([1,3]/10)
psum = p1 + p2 - p3
@test degree(psum) == 1         # will have wrong degree
@test degree(truncate(psum)) == 0 # the degree should be correct after truncation

@test truncate(Poly([2,1]),rtol=1/2,atol=0) == Poly([2])
@test truncate(Poly([2,1]),rtol=1,atol=0)   == Poly([0])
@test truncate(Poly([2,1]),rtol=0,atol=1)   == Poly([2])

@test norm(Poly([1., 2.])) == norm([1., 2.])
@test norm(Poly([1., 2.]), 1) == norm([1., 2.], 1)

## setindex!
println("Test for setindex!()")
p1    = Poly([1,2,1])
p1[5] = 1
@test p1[5] == 1
@test p1 == Poly([1,2,1,0,0,1])

p1[:] = 0
@test p1 ≈ zero(p1)
@test zero(Poly{Int}) == Poly(Int[])
@test one(Poly{Int}) == Poly([1])
p1[:] = [1,2,1,0,0,1]
@test p1 == Poly([1,2,1,0,0,1])

## elementwise operations #52
println("Test for element-wise operations")
p1  = Poly([1, 2])
p2  = Poly([3, 1.])
p   = [p1, p2]
q   = [3, p1]
@test isa(q,Vector{Poly{Int64}})
psum  = p .+ 3
pprod = p .* 3
pmin  = p .- 3
@test isa(psum, Vector{Poly{Float64}})
@test isa(pprod,Vector{Poly{Float64}})
@test isa(pmin, Vector{Poly{Float64}})

## getindex with ranges #43
p1 = Poly([4,5,6])
@test all(p1[0:1] .== [4,5])
@test all(p1[0:end] .== [4,5,6])
p1[0:1] = [7,8]
@test all(p1[0:end] .== [7,8,6])

@test p1[:] == coeffs(p1)

## conjugate of poly (issue #59)
as = [im, 1, 2]
bs = [1, 1, 2]
@test conj(Poly(as)) == Poly(conj(as))
@test conj(Poly(bs)) == Poly(conj(bs))
## and transpose
@test transpose(Poly(as)) == Poly(as)

## unnecessary copy in convert #65
p1 = Poly([1,2])
p2 = convert(Poly{Int64}, p1)
p2[3] = 3
@test p1[3] == 3

## Polynomials with non-Real type
import Base: +, *, -
struct Mod2 <: Number
  v::Bool
end
+(x::Mod2,y::Mod2) = Mod2(xor(x.v, y.v))
*(x::Mod2,y::Mod2) = Mod2(x.v&y.v)
-(x::Mod2,y::Mod2) = x+y
-(x::Mod2) = x
Base.one(::Type{Mod2}) = Mod2(true)
Base.zero(::Type{Mod2}) = Mod2(false)
Base.convert(::Type{Mod2},x::Integer) = Mod2(convert(Bool,x))
Base.convert(::Type{Bool},x::Mod2) = x.v

# Test that none of this throws
p = Poly([Mod2(true),Mod2(false), Mod2(true)])
repr(p)
@test p(Mod2(false)) == p[0]


## changes to show
p = Poly([1,2,3,1])  # leading coefficient of 1
@test repr(p) == "Poly(1 + 2*x + 3*x^2 + x^3)"
p = Poly([1.0, 2.0, 3.0, 1.0])
@test repr(p) == "Poly(1.0 + 2.0*x + 3.0*x^2 + 1.0*x^3)"
p = Poly([1, im])
@test repr(p) == "Poly(1 + im*x)"
p = Poly([1+im, 1-im, -1+im, -1 - im])# minus signs
@test repr(p) == "Poly((1 + 1im) + (1 - 1im)*x - (1 - 1im)*x^2 - (1 + 1im)*x^3)"
p = Poly([1.0, 0 + NaN*im, NaN, Inf, 0 - Inf*im]) # handle NaN or Inf appropriately
@test repr(p) == "Poly(1.0 + NaN*im*x + NaN*x^2 + Inf*x^3 - Inf*im*x^4)"

p = Poly([1,2,3])

if VERSION < v"0.7-"
    @test reprmime("text/latex", p) == "\$1 + 2\\cdot x + 3\\cdot x^{2}\$"
    p = Poly([1//2, 2//3, 1])
    @test reprmime("text/latex", p) == "\$\\frac{1}{2} + \\frac{2}{3}\\cdot x + x^{2}\$"
else
    @test repr("text/latex", p) == "\$1 + 2\\cdot x + 3\\cdot x^{2}\$"
    p = Poly([1//2, 2//3, 1])
    @test repr("text/latex", p) == "\$\\frac{1}{2} + \\frac{2}{3}\\cdot x + x^{2}\$"
end

# customized printing with printpoly
function printpoly_to_string(args...; kwargs...)
    buf = IOBuffer()
    printpoly(buf, args...; kwargs...)
    return Compat.String(take!(buf))
end
@test printpoly_to_string(Poly([1,2,3], "y")) == "1 + 2*y + 3*y^2"
@test printpoly_to_string(Poly([1,2,3], "y"), descending_powers=true) == "3*y^2 + 2*y + 1"

## want to be able to copy and paste
if VERSION < v"0.6.0-dev"
    string_eval_poly(p,x) = eval(Expr(:function, Expr(:call, :f, :x), parse(string(p)[6:end-1])))(x)
    p = Poly([1,2,3]) # copy and paste
    q = Poly([1//1, 2//1, 3//1])
    r = Poly([1.0, 2, 3])
    @test string_eval_poly(p, 5) == p(5)
    @test string_eval_poly(q, 5) == q(5)
    @test string_eval_poly(r, 5) == r(5)
end
## check hashing
p = poly([1,2,3])
q = poly([1,2,3])
@test hash(p) == hash(q)

## Check for Inf/NaN operations
p1 = Poly([Inf, Inf])
p2 = Poly([0, Inf])
@test p1(Inf) == Inf
@test isnan(p1(-Inf))
@test isnan(p1(0))
@test p2(-Inf) == -Inf

## Check for isequal
p1 = Poly([-0., 5., Inf])
p2 = Poly([0., 5., Inf])
p3 = Poly([0, NaN])

@test p1 == p2 && !isequal(p1, p2)
@test p3 === p3 && p3 ≠ p3 && isequal(p3, p3)

## Handling of `NaN`s
p     = Poly([NaN, 1, 5])
pder  = polyder(p)
pint  = polyint(p)

@test isnan(p(1))                 # p(1) evaluates to NaN
@test isequal(pder, Poly([NaN]))  # derivative will give Poly([NaN])
@test isequal(pint, Poly([NaN]))  # integral will give Poly([NaN])

pint  = polyint(p, Complex(0.))
@test isequal(pint, Poly([NaN]))  # integral will give Poly([NaN])

## proper conversions in arithmetic with different element-types #94
p = Poly([0,one(Float64)])
@test Poly{Complex{Float64}} == typeof(p+1im)
@test Poly{Complex{Float64}} == typeof(1im-p)
@test Poly{Complex{Float64}} == typeof(p*1im)

## comparison between `Number`s and `Poly`s
p1s = Poly([1,2], :s)
p1x = Poly([1,2], :x)
p2s = Poly([1], :s)

@test p1s == p1s
@test p1s ≠ p1x
@test p1s ≠ p2s

@test_throws ErrorException p1s ≈ p1x
@test p1s ≉ p2s
@test p1s ≈ Poly([1,2.], :s)

@test p2s ≈ 1. ≈ p2s
@test p2s == 1. == p2s
@test p2s ≠ 2. ≠ p2s
@test p1s ≠ 2. ≠ p1s

@test nnz(map(Poly, sparse(1.0I, 5,5))) == 5

@test Poly([0.5]) + 2 == Poly([2.5])
@test 2 - Poly([0.5]) == Poly([1.5])

# test size
@test size(Poly([0.5, 0.2])) == (2,)
@test size(Poly([0.5, 0.2]), 1) == 2
@test size(Poly([0.5, 0.2]), 1, 2) == (2,1)

# test iteration
p1 = Poly([1,2,0,3])
for term in p1
  @test isa(term, Poly)
end

@test eltype(typeof(p1)) == typeof(p1)

# after implementing the iteration interface, i.e., `start` and its friends,
# previously throwing `collect` function (and maybe, other similar functions)
# started giving different errors due to `eltype{T}(::Poly{T}) = T`. Changing
# this definition has resulted in JuliaLang/METADATA.jl#8528. One small fix
# was to direct `collect{T}(p::Poly{T})` to `collect(Poly{T}, p)`.

@test eltype(p1) == Int
@test eltype(collect(p1)) == Poly{Int}
@test eltype(collect(Poly{Float64}, p1)) == Poly{Float64}
@test_throws InexactError collect(Poly{Int}, Poly([1.2]))

@test length(collect(p1)) == degree(p1)+1

@test [p1[idx] for idx in eachindex(p1)] == [1,2,0,3]

p1 = Poly([2.,5.,1.])
p2 = Poly([1.,2.,3.])

@test degree(gcd(p1, p2))                       == 0          # no common roots
@test degree(gcd(p1, Poly(5)))                  == 0          # ditto
@test degree(gcd(p1, Poly(eps(0.))))            == 0          # ditto
@test degree(gcd(p1, Poly(0)))                  == degree(p1) # Poly(0) has the roots of p1
@test degree(gcd(p1+p2*170.10734737144486, p2)) == 0          # see, c.f., #122

p1 = poly([1.,2.,3.])
p2 = poly([1.,2.,6.])

@test (res = roots(gcd(p1, p2)); 1. ∈ res && 2. ∈ res)
