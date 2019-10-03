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
    @test typeof(p).parameters[1] == eltype(coeff)
    @test eltype(p) == eltype(coeff)
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
    p0 = Polynomial{Float64}([0])
    p1 = Polynomial{Float64}([1])
    p2 = Polynomial{Float64}([5, 6, -3, 2 ,4])
    p3 = Polynomial{Float64}([7, -3, 2, 6])
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

@testset "Identities" begin
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
end

@testset "Fitting" begin
    xs = range(0, stop = pi, length = 10)
    ys = sin.(xs)

    p = fit(xs, ys)
    y_fit = p.(xs)
    abs_error = abs.(y_fit .- ys)
    @test maximum(abs_error) <= 0.03

    p = fit(xs, ys, 2)
    y_fit = p.(xs)
    abs_error = abs.(y_fit .- ys)
    @test maximum(abs_error) <= 0.03
end

@testset "Values" begin
    @test fromroots(Int[])(2.) == 1.
    @test pN(-.125) == 276.9609375
    @test pNULL(10) == 0
    @test p0(-10) == 0
    @test fromroots([1 // 2, 3 // 2])(1 // 2) == 0 // 1
end

@testset "Conversion" begin
    
end

@testset "Roots" begin
    # From roots
    r = [3, 2]
    p = fromroots(Polynomial, r)
    @test fromroots(r) == Polynomial([6, -5, 1])
    @test p == Polynomial([6, -5, 1])
    @test roots(p) ≈ r

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
    int = integral(der, 1)
    @test int.coeffs == c

    @test_broken derivative(integral(pN)) == pN
    @test derivative(pR) == Polynomial([-2 // 1,2 // 1])
    @test derivative(p3) == Polynomial([2,2])
    @test derivative(p1) == derivative(p0) == derivative(pNULL) == pNULL
    @test_throws ErrorException derivative(pR, -1)
    @test_broken derivative(pNULL, 1) == p1
    @test_broken derivative(Polynomial(Rational[1,2,3])) == Polynomial(Rational[0, 1, 1, 1])
    @test integrate(Polynomial([1,1,0,0]), 0, 2) == 4.0

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
end

@testset "Linear Algebra" begin
    p = Polynomial([3, 4])
    @test norm(p) == 5
    p = Polynomial([-1, 3, 5, -2])
    @test norm(p) ≈ 6.244997998398398
    p = Polynomial([1 - 1im, 2 - 3im])
    p2 = conj(p)
    @test p2.coeffs == [1 + 1im, 2 + 3im]
end

@testset "Indexing" begin
    # Indexing
    p = Polynomial([-1, 3, 5, -2])
    @test p[0] == -1
    @test p[[1, 2]] == [3, 5]
    @test p[1:2] == [3, 5]
    @test p[:] == [-1, 3, 5, -2]

    @test p[end] == p.coeffs[end]
    p[1] = 2
    @test p.coeffs[2] == 2
    p[2:3] = [1, 2]
    @test p.coeffs[3:4] == [1, 2]
    p[0:1] = 0
    @test p.coeffs[1:2] == [0, 0]
    p[:] = 1
    @test p.coeffs == ones(4)
    
end

