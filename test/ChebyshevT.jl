@testset "Construction" for coeff in [
    Int64[1, 1, 1, 1],
    Float32[1, -4, 2],
    ComplexF64[1 - 1im, 2 + 3im],
    [3 // 4, -2 // 1, 1 // 1]
]  
    p = ChebyshevT(coeff)
    @test p.coeffs == coeff
    @test coeffs(p) == coeff
    @test degree(p) == length(coeff) - 1
    @test order(p) == length(p) == length(coeff)
    @test p.var == :x
    @test length(p) == length(coeff)
    @test size(p) == size(coeff)
    @test size(p, 1) == size(coeff, 1)
    @test typeof(p).parameters[1] == eltype(coeff)
    @test eltype(p) == eltype(coeff)
end

@testset "Other Construction" begin
    # Leading 0s
    p = ChebyshevT([1, 2, 0, 0])
    @test p.coeffs == [1, 2]
    @test length(p) == 2

    # Different type
    p = ChebyshevT{Float64}(ones(Int32, 4))
    @test p.coeffs == ones(Float64, 4)

    p = ChebyshevT(30)
    @test p.coeffs == [30]

    p = zero(ChebyshevT{Int})
    @test p.coeffs == [0]
    
    p = one(ChebyshevT{Int})
    @test p.coeffs == [1]

    pNULL = ChebyshevT(Int[])
    @test iszero(pNULL)
    @test degree(pNULL) == -1

    p0 = ChebyshevT([0])
    @test iszero(p0)
    @test degree(p0) == -1
end


@testset "Values" begin
    c1 = ChebyshevT([2.5, 2.0, 1.5]) #1 + 2x + 3x^2
    x = -1:0.5:1
    expected = [2.0, 0.75, 1.0, 2.75, 6.0]
    @test c1.(x) == expected

    @test c1 == Polynomial([1, 2, 3])
end
