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
end

@testset "Conversions" begin
    c1 = ChebyshevT([2.5, 2.0, 1.5])
    @test c1 ≈ Polynomial([1, 2, 3])
    @test convert(Polynomial{Float64}, c1) == Polynomial{Float64}([1, 2, 3])

end

@testset "Arithmetic" begin
    # multiplication
    for i in 1:5, j in 1:5
        target = zeros(i + j + 1)
        target[end] += 0.5
        target[abs(i - j) + 1] += 0.5
        c1 = ChebyshevT(append!(zeros(i), [1]))
        c2 = ChebyshevT(append!(zeros(j), [1]))
        @test c1 * c2 ≈ ChebyshevT(target)
    end

    c1 = ChebyshevT([1, 2, 3])
    c2 = ChebyshevT([3, 2, 1])
    @test c1 * c2 == ChebyshevT([6.5, 12, 12, 4, 1.5])

end

@testset "z-series" for i in 1:5
    # c to z
    input = append!([2], ones(i))
    target = append!(append!(0.5 .* ones(i), 2), 0.5 .* ones(i))
    zs = Polynomials._c_to_z(input)
    @test zs == target
    c = Polynomials._z_to_c(zs)
    @test c == input
end
