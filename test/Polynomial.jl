using LinearAlgebra

@testset "Construction" for coeff in [
    Int64[1, 1, 1, 1],
    Float32[1, -4, 2],
    ComplexF64[1 - 1im, 2 + 3im],
    [3 // 4, -2 // 1, 1 // 1]
]
    p = Polynomial(coeff)
    @test p.coeffs == coeff
    @test coeffs(p) == coeff
    @test degree(p) == length(coeff) - 1
    @test p.var == :x
    @test length(p) == length(coeff)
    @test size(p) == size(coeff)
    @test size(p, 1) == size(coeff, 1)
    @test typeof(p).parameters[1] == eltype(coeff)
    @test eltype(p) == eltype(coeff)
    @test all([-200, -0.3, 1, 48.2] .∈ domain(p))
end

@testset "Mapdomain" begin
    x = -30:20
    mx = mapdomain(Polynomial, x)
    @test mx == x

    x = 0.5:0.01:0.6
    mx = mapdomain(Polynomial, x)
    @test mx == x
end

@testset "Other Construction" begin
    # Leading 0s
    p = Polynomial([1, 2, 0, 0])
    @test p.coeffs == [1, 2]
    @test length(p) == 2

    # Different type
    p = Polynomial{Float64}(ones(Int32, 4))
    @test p.coeffs == ones(Float64, 4)

    p = Polynomial(30)
    @test p.coeffs == [30]

    p = zero(Polynomial{Int})
    @test p.coeffs == [0]

    p = one(Polynomial{Int})
    @test p.coeffs == [1]

    pNULL = Polynomial(Int[])
    @test iszero(pNULL)
    @test degree(pNULL) == -1

    p0 = Polynomial([0])
    @test iszero(p0)
    @test degree(p0) == -1
end

pNULL = Polynomial(Int[])
p0 = Polynomial([0])
p1 = Polynomial([1,0,0,0,0,0,0,0,0,0,0,0,0,0])
p2 = Polynomial([1,1,0,0])
p3 = Polynomial([1,2,1,0,0,0,0])
p4 = Polynomial([1,3,3,1,0,0])
p5 = Polynomial([1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
pN = Polynomial([276,3,87,15,24,0])
pR = Polynomial([3 // 4, -2 // 1, 1 // 1])

@testset "Arithmetic" begin
    @test p3 == Polynomial([1,2,1])
    @test pN * 10 == Polynomial([2760, 30, 870, 150, 240])
    @test pN / 10.0 == Polynomial([27.6, 0.3, 8.7, 1.5, 2.4])
    @test 10 * pNULL + pN == pN
    @test 10 * p0 + pN == pN
    @test p5 + 2 * p1 == Polynomial([3,4,6,4,1])
    @test 10 * pNULL - pN == -pN
    @test p0 - pN == -pN
    @test p5 - 2 * p1 == Polynomial([-1,4,6,4,1])
    @test p2 * p2 * p2 == p4
    @test p2^4 == p5
    @test pNULL^3 == pNULL
    @test pNULL * pNULL == pNULL

    @test pNULL + 2 == p0 + 2 == 2 + p0 == Polynomial([2])
    @test p2 - 2 == -2 + p2 == Polynomial([-1,1])
    @test 2 - p2 == Polynomial([1,-1])
end

@testset "Divrem" begin
    p0 = Polynomial([0])
    p1 = Polynomial([1])
    p2 = Polynomial([5, 6, -3, 2 ,4])
    p3 = Polynomial([7, -3, 2, 6])
    p4 = p2 * p3
    @test divrem(p4, p2) == (p3, zero(p3))
    @test p3 % p2 == p3
    @test all((map(abs, (p2 ÷ p3 - Polynomial([1 / 9,2 / 3])).coeffs)) .< eps())
    @test divrem(p0, p1) == (p0, p0)
    @test divrem(p1, p1) == (p1, p0)
    @test divrem(p2, p2) == (p1, p0)
    @test divrem(pR, pR) == (one(pR), zero(pR))
    @test_throws DivideError p1 ÷ p0
    @test_throws DivideError divrem(p0, p0)
end

@testset "Comparisons" begin
    pX = Polynomial([1, 2, 3, 4, 5])
    pS1 = Polynomial([1, 2, 3, 4, 5], "s")
    pS2 = Polynomial([1, 2, 3, 4, 5], 's')
    pS3 = Polynomial([1, 2, 3, 4, 5], :s)
    @test pX != pS1
    @test pS1 == pS2
    @test pS1 == pS3
    @test_throws ErrorException pS1 + pX
    @test_throws ErrorException pS1 - pX
    @test_throws ErrorException pS1 * pX
    @test_throws ErrorException pS1 ÷ pX
    @test_throws ErrorException pS1 % pX

    # Testing copying.
    pcpy1 = Polynomial([1,2,3,4,5], :y)
    pcpy2 = copy(pcpy1)
    @test pcpy1 == pcpy2

    # Check for isequal
    p1 = Polynomial([-0., 5., Inf])
    p2 = Polynomial([0., 5., Inf])
    p3 = Polynomial([0, NaN])

    @test p1 == p2 && !isequal(p1, p2)
    @test p3 === p3 && p3 ≠ p3 && isequal(p3, p3)

    p = fromroots(Polynomial, [1,2,3])
    q = fromroots(Polynomial, [1,2,3])
    @test hash(p) == hash(q)

    p1s = Polynomial([1,2], :s)
    p1x = Polynomial([1,2], :x)
    p2s = Polynomial([1], :s)

    @test p1s == p1s
    @test p1s ≠ p1x
    @test p1s ≠ p2s

    @test_throws ErrorException p1s ≈ p1x
    @test p1s ≉ p2s
    @test p1s ≈ Polynomial([1,2.], :s)

    @test p2s ≈ 1.0 ≈ p2s
    @test p2s == 1.0 == p2s
    @test p2s ≠ 2.0 ≠ p2s
    @test p1s ≠ 2.0 ≠ p1s

    @test nnz(map(Polynomial, sparse(1.0I, 5, 5))) == 5

    @test Polynomial([0.5]) + 2 == Polynomial([2.5])
    @test 2 - Polynomial([0.5]) == Polynomial([1.5])
end

@testset "Fitting" begin
    xs = range(0, stop = π, length = 10)
    ys = sin.(xs)

    p = fit(Polynomial, xs, ys)
    y_fit = p.(xs)
    abs_error = abs.(y_fit .- ys)
    @test maximum(abs_error) <= 0.03

    p = fit(Polynomial, xs, ys, 2)
    y_fit = p.(xs)
    abs_error = abs.(y_fit .- ys)
    @test maximum(abs_error) <= 0.03

    # Test weighted
    for W in [1, ones(size(xs)), diagm(0 => ones(size(xs)))]
        p = fit(Polynomial, xs, ys, 2, weights = W)
        @test p.(xs) ≈ y_fit
    end


    # Getting error on passing Real arrays to polyfit #146
    xx = Real[20.0, 30.0, 40.0]
    yy = Real[15.7696, 21.4851, 28.2463]
    fit(xx, yy, 2)
end

@testset "Values" begin
    @test fromroots(Polynomial, Int[])(2.) == 1.
    @test pN(-.125) == 276.9609375
    @test pNULL(10) == 0
    @test p0(-10) == 0
    @test fromroots(Polynomial, [1 // 2, 3 // 2])(1 // 2) == 0 // 1

    # Check for Inf/NaN operations
    p1 = Polynomial([Inf, Inf])
    p2 = Polynomial([0, Inf])
    @test p1(Inf) == Inf
    @test isnan(p1(-Inf))
    @test isnan(p1(0))
    @test p2(-Inf) == -Inf

    # issue #189
    p = Polynomial([0,1,2,3])
    A = [0 1; 0  0];
    @test  p(A) == A  + 2A^2 + 3A^3

end

@testset "Conversion" begin
    # unnecessary copy in convert #65
    p1 = Polynomial([1,2])
    p2 = convert(Polynomial{Int64}, p1)
    p2[3] = 3
    @test p1[3] == 3

    p = Polynomial([0,one(Float64)])
    @test Polynomial{Complex{Float64}} == typeof(p + 1im)
    @test Polynomial{Complex{Float64}} == typeof(1im - p)
    @test Polynomial{Complex{Float64}} == typeof(p * 1im)
end

@testset "Roots" begin
    # From roots
    r = [2, 3]
    p = fromroots(Polynomial, r)
    @test fromroots(r) == Polynomial([6, -5, 1])
    @test p == Polynomial([6, -5, 1])
    @test sort(roots(p)) ≈ r

    @test roots(p0) == roots(p1) == roots(pNULL) == []
    @test roots(Polynomial([0,1,0])) == [0.0]
    r = roots(Polynomial([0,1,0]))

    @test roots(p2) == [-1]
    a_roots = copy(pN.coeffs)
    @test all(map(abs, sort(roots(fromroots(a_roots))) - sort(a_roots)) .< 1e6)
    @test length(roots(p5)) == 4
    @test roots(pNULL) == []
    @test sort(roots(pR)) == [1 // 2, 3 // 2]

    A = [1 0; 0 1]
    p = fromroots(Polynomial, A)
    @test fromroots(A) == Polynomial(Float64[1, -2, 1])
    @test p == Polynomial(Float64[1, -2, 1])
    @test roots(p) ≈ sort!(eigvals(A), rev = true)

    x = variable()
    plarge = 8.362779449448982e41 - 2.510840694154672e57x + 4.2817430781178795e44x^2 - 1.6225927682921337e31x^3 + 1.0x^4  # #120
    @test length(roots(plarge)) == 4
end

@testset "Integrals and Derivatives" begin
    # Integrals derivatives
    c = [1, 2, 3, 4]
    p = Polynomial(c)
    der = derivative(p)
    @test der.coeffs == [2, 6, 12]
    int = integrate(der, 1)
    @test int.coeffs == c

    @test derivative(integrate(pN)) == convert(Polynomial{Float64}, pN)
    @test derivative(pR) == Polynomial([-2 // 1,2 // 1])
    @test derivative(p3) == Polynomial([2,2])
    @test derivative(p1) == derivative(p0) == derivative(pNULL) == pNULL
    @test_throws ErrorException derivative(pR, -1)
    @test integrate(pNULL, 1) == convert(Polynomial{Float64}, p1)
    rc = Rational{Int64}[1,2,3]
    @test integrate(Polynomial(rc)) == Polynomial{eltype(rc)}([0, 1, 1, 1])
    @test integrate(Polynomial([1,1,0,0]), 0, 2) == 4.0

    for i in 1:10
        p = Polynomial{Float64}(rand(1:5, 6))
        @test degree(round(p - integrate(derivative(p)), digits=13)) <= 0
        @test degree(round(p - derivative(integrate(p)), digits=13)) <= 0
    end


    # Handling of `NaN`s
    p     = Polynomial([NaN, 1, 5])
    pder  = derivative(p)
    pint  = integrate(p)

    @test isnan(p(1)) # p(1) evaluates to NaN
    @test isequal(pder, Polynomial([NaN]))
    @test isequal(pint, Polynomial([NaN]))

    pint  = integrate(p, 0.0im)
    @test isequal(pint, Polynomial([NaN]))

    # Issue with overflow and polyder Issue #159
    @test !iszero(derivative(Polynomial(BigInt[0, 1])^100, 100))
end

@testset "Elementwise Operations" begin
    p1  = Polynomial([1, 2])
    p2  = Polynomial([3, 1.])
    p   = [p1, p2]
    q   = [3, p1]
    @test q isa Vector{Polynomial{Int64}}
    psum  = p .+ 3
    pprod = p .* 3
    pmin  = p .- 3
    @test psum  isa Vector{Polynomial{Float64}}
    @test pprod isa Vector{Polynomial{Float64}}
    @test pmin  isa Vector{Polynomial{Float64}}
end

@testset "Chop and Truncate" begin
    # chop and truncate
    p = Polynomial([1, 1, 1, 1])
    p.coeffs[end] = 0
    @assert p.coeffs == [1, 1, 1, 0]
    p = chop(p)
    @test p.coeffs == [1, 1, 1]
    ## truncation
    p1 = Polynomial([1,1] / 10)
    p2 = Polynomial([1,2] / 10)
    p3 = Polynomial([1,3] / 10)
    psum = p1 + p2 - p3
    @test degree(psum) == 1         # will have wrong degree
    @test degree(truncate(psum)) == 0 # the degree should be correct after truncation

    @test truncate(Polynomial([2,1]), rtol = 1 / 2, atol = 0) == Polynomial([2])
    @test truncate(Polynomial([2,1]), rtol = 1, atol = 0)   == Polynomial([0])
    @test truncate(Polynomial([2,1]), rtol = 0, atol = 1)   == Polynomial([2])

    pchop = Polynomial([1, 2, 3, 0, 0, 0])
    pchopped = chop(pchop)
    @test roots(pchop) == roots(pchopped)

    # round
    psmall = Polynomial(eps()*rand(1:10,  9))
    @test degree(round(psmall, digits=14)) == -1

end

@testset "Linear Algebra" begin
    p = Polynomial([3, 4])
    @test norm(p) == 5
    p = Polynomial([-1, 3, 5, -2])
    @test norm(p) ≈ 6.244997998398398
    p = Polynomial([1 - 1im, 2 - 3im])
    p2 = conj(p)
    @test p2.coeffs == [1 + 1im, 2 + 3im]
    @test transpose(p) == p
    @test transpose!(p) == p

    @test norm(Polynomial([1., 2.])) == norm([1., 2.])
    @test norm(Polynomial([1., 2.]), 1) == norm([1., 2.], 1)
end

@testset "Indexing" begin
    # Indexing
    p = Polynomial([-1, 3, 5, -2])
    @test p[0] == -1
    @test p[[1, 2]] == [3, 5]
    @test p[1:2] == [3, 5]
    @test p[:] == [-1, 3, 5, -2]

    p1    = Poly([1,2,1])
    p1[5] = 1
    @test p1[5] == 1
    @test p1 == Polynomial([1,2,1,0,0,1])

    @test p[end] == p.coeffs[end]
    p[1] = 2
    @test p.coeffs[2] == 2
    p[2:3] = [1, 2]
    @test p.coeffs[3:4] == [1, 2]
    p[0:1] = 0
    @test p.coeffs[1:2] == [0, 0]
    p[:] = 1
    @test p.coeffs == ones(4)

    p[:] = 0
    @test chop(p) ≈ zero(p)

    p1 = Polynomial([1,2,0,3])
    for term in p1
        @test isa(term, Polynomial)
    end

    @test eltype(p1) == Int
    @test eltype(collect(p1)) == Polynomial{Int}
    @test eltype(collect(Polynomial{Float64}, p1)) == Polynomial{Float64}
    @test_throws InexactError collect(Polynomial{Int}, Polynomial([1.2]))

    @test length(collect(p1)) == degree(p1) + 1

    @test [p1[idx] for idx in eachindex(p1)] == [1,2,0,3]
end

@testset "Copying" begin
    pcpy1 = Polynomial([1,2,3,4,5], :y)
    pcpy2 = copy(pcpy1)
    @test pcpy1 == pcpy2
end

@testset "GCD" begin
    p1 = Polynomial([2.,5.,1.])
    p2 = Polynomial([1.,2.,3.])

    @test degree(gcd(p1, p2)) == 0          # no common roots
    @test degree(gcd(p1, Polynomial(5))) == 0          # ditto
    @test degree(gcd(p1, Polynomial(eps(0.)))) == 0          # ditto
    @test degree(gcd(p1, Polynomial(0))) == degree(p1) # Polynomial(0) has the roots of p1
    @test degree(gcd(p1 + p2 * 170.10734737144486, p2)) == 0          # see, c.f., #122

    p1 = fromroots(Polynomial, [1.,2.,3.])
    p2 = fromroots(Polynomial, [1.,2.,6.])
    res = roots(gcd(p1, p2))
    @test 1. ∈ res
    @test 2. ∈ res
end

@testset "Pade" begin
    # exponential
    coeffs = 1 .// BigInt.(gamma.(1:17))
    a = Polynomial(coeffs)
    PQexp = Pade(a, 8, 8)
    @test PQexp(1.0) ≈ exp(1.0)
    @test PQexp(-1.0) ≈ exp(-1.0)

    # sine
    coeffs = BigInt.(sinpi.((0:16) ./ 2)) .// BigInt.(gamma.(1:17))
    p = Polynomial(coeffs)
    PQsin = Pade(p, 8, 7)
    @test PQsin(1.0) ≈ sin(1.0)
    @test PQsin(-1.0) ≈ sin(-1.0)

    # cosine
    coeffs = BigInt.(sinpi.((1:17) ./ 2)) .// BigInt.(gamma.(1:17))
    p = Polynomial(coeffs)
    PQcos = Pade(p, 8, 8)
    @test PQcos(1.0) ≈ cos(1.0)
    @test PQcos(-1.0) ≈ cos(-1.0)

    # summation of a factorially divergent series
    γ = 0.5772156649015
    s = BigInt.(gamma.(BigInt(1):BigInt(61)))
    coeffs = (BigInt(-1)).^(0:60) .* s .// 1
    d = Polynomial(coeffs)
    PQexpint = Pade(d, 30, 30)
    @test Float64(PQexpint(1.0)) ≈ exp(1) * (-γ - sum([(-1)^k / k / gamma(k + 1) for k = 1:20]))
end

@testset "Showing" begin
    p = Polynomial([1, 2, 3])
    @test sprint(show, p) == "Polynomial(1 + 2*x + 3*x^2)"

    p = Polynomial([1.0, 2.0, 3.0])
    @test sprint(show, p) == "Polynomial(1.0 + 2.0*x + 3.0*x^2)"

    p = Polynomial([1 + 1im, -2im])
    @test sprint(show, p) == "Polynomial((1 + 1im) - 2im*x)"

    p = Polynomial{Rational}([1, 4])
    @test sprint(show, p) == "Polynomial(1//1 + 4//1*x)"

    p = Polynomial([1,2,3,1])  # leading coefficient of 1
    @test repr(p) == "Polynomial(1 + 2*x + 3*x^2 + x^3)"
    p = Polynomial([1.0, 2.0, 3.0, 1.0])
    @test repr(p) == "Polynomial(1.0 + 2.0*x + 3.0*x^2 + 1.0*x^3)"
    p = Polynomial([1, im])
    @test repr(p) == "Polynomial(1 + im*x)"
    p = Polynomial([1 + im, 1 - im, -1 + im, -1 - im])# minus signs
    @test repr(p) == "Polynomial((1 + 1im) + (1 - 1im)*x - (1 - 1im)*x^2 - (1 + 1im)*x^3)"
    p = Polynomial([1.0, 0 + NaN * im, NaN, Inf, 0 - Inf * im]) # handle NaN or Inf appropriately
    @test repr(p) == "Polynomial(1.0 + NaN*im*x + NaN*x^2 + Inf*x^3 - Inf*im*x^4)"

    p = Polynomial([1,2,3])

    @test repr("text/latex", p) == "\$1 + 2\\cdot x + 3\\cdot x^{2}\$"
    p = Polynomial([1 // 2, 2 // 3, 1])
    @test repr("text/latex", p) == "\$\\frac{1}{2} + \\frac{2}{3}\\cdot x + x^{2}\$"

    # customized printing with printpoly
    function printpoly_to_string(args...; kwargs...)
        buf = IOBuffer()
        printpoly(buf, args...; kwargs...)
        return String(take!(buf))
    end
    @test printpoly_to_string(Polynomial([1,2,3], "y")) == "1 + 2*y + 3*y^2"
    @test printpoly_to_string(Polynomial([1,2,3], "y"), descending_powers = true) == "3*y^2 + 2*y + 1"
    @test printpoly_to_string(Polynomial([2, 3, 1], :z), descending_powers = true, offset = -2) == "1 + 3*z^-1 + 2*z^-2"
    @test printpoly_to_string(Polynomial([-1, 0, 1], :z), offset = -1, descending_powers = true) == "z - z^-1"
end

@testset "Plotting" begin
    p = fromroots([-1, 1]) # x^2 - 1
    r = -1.4:0.055999999999999994:1.4
    rec = apply_recipe(Dict{Symbol,Any}(), p)
    @test getfield(rec[1], 1) == Dict{Symbol,Any}(:label => "-1 + x^2")
    @test rec[1].args == (r, p.(r))

    r = -1:0.02:1
    rec = apply_recipe(Dict{Symbol,Any}(), p, -1, 1)
    @test getfield(rec[1], 1) == Dict{Symbol,Any}(:label => "-1 + x^2")
    @test rec[1].args == (r, p.(r))

    p = ChebyshevT([1,1,1])
    rec = apply_recipe(Dict{Symbol,Any}(), p)
    r = -1.0:0.02:1.0  # uses domain(p)
    @test rec[1].args == (r, p.(r))

end
