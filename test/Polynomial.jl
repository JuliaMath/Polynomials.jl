@testset "Construction" for coeff in [
    Int64[1, 1, 1, 1],
    Float32[1, -4, 2],
    ComplexF64[1 - 1im, 2 + 3im]
]  
    p = Polynomial(coeff)
    @test p.coeffs == coeff
    @test p.var == :x
    @test length(p) == length(coeff)
    @test typeof(p).parameters[1] == eltype(coeff)
end

@testset "From Roots" begin
    r = [3, 2]
    @test roots(r, Polynomial) == Polynomial([6, -5, 1])
    A = [1 0; 0 1]
    @test roots(A, Polynomial) == Polynomial(Float64[1, -2, 1])
end

