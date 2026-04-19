using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(Polynomials)
end

if VERSION ≥ v"1.12.0"
    @testset "JET" begin
        JET.test_package(Polynomials, ignored_modules=(AnyFrameModule(LinearAlgebra), AnyFrameModule(Base)))
    end
end
