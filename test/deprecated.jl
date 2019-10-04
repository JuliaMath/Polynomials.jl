
@testset "Construction" begin
    @test_deprecated Poly(1)
    @test_deprecated Poly(2, "x")
    @test_deprecated Poly(collect(1:10))
    @test_deprecated Poly(collect(1:10), "x")
    @test_deprecated Poly{Float64}(collect(1:10))
    @test_deprecated Poly{Float64}(collect(1:10), "x")
    p = Polynomial([1, 2, 3])
    @test_deprecated p.a
end

@testset "From roots" begin
    @test_deprecated poly([1, 2])
    @test_deprecated poly([1, 2], "s")

    M = [1 2; 3 4]
    @test_deprecated poly(M)
    @test_deprecated poly(M, "s")
end

@testset "Evaluation" begin
    p1 = Polynomial([1, 2, 3, 0, 5])

    @test_deprecated polyval(p1, 1)
    @test_deprecated polyval(p1, randn(10))
    @test_deprecated p1(randn(10))
end

@testset "Integrals and Derivatives" begin
    p1 = Polynomial([1, 2, 3, 0, 5])
    @test_deprecated polyint(p1)
    @test_deprecated polyint(p1, 3)
    @test_deprecated polyint(p1, 0, 5)

    @test_deprecated polyder(p1)
    @test_deprecated polyder(p1, 2)
end

@testset "Fits and Pade" begin
    xs = range(0, stop = π, length = 10)
    ys = sin.(xs)

    @test_deprecated polyfit(xs, ys)
    @test_deprecated polyfit(xs, ys, 2)
    @test_deprecated polyfit(xs, ys, :t)

    a = Polynomial(1:17)
    PQexp = Pade(a, 8, 8)
    @test_deprecated padeval(PQexp, 1.0)
end

@testset "Converts to Polynomial" begin
    p1 = Poly([1, 2, 3])
    p = convert(Polynomial, p1)
    @test p.coeffs == p1.coeffs
    @test p.var == p1.var

    p = convert(Polynomial{Float64}, p1)
    @test p.coeffs == Float64[1, 2, 3]
    @test p.var == p1.var

    @test p1 + p == Polynomial([2, 4, 6])
end