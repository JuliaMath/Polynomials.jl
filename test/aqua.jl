using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(Polynomials)
end

# Any two-argument `Base.promote` method, whatever its signature, invalidates all
# compiled `promote(a, b)` call sites whose arguments are inferred abstractly.
# Promotion should go through `promote_rule` and `convert` instead.
@testset "No two-argument promote methods" begin
    ms = [m for m in methods(promote)
              if Base.moduleroot(m.module) === Polynomials && m.nargs == 3 && !m.isva]
    @test isempty(ms)
end

# JET is not functional on unreleased Julia versions (e.g. nightly)
if VERSION ≥ v"1.12.0" && isempty(VERSION.prerelease)
    @testset "JET" begin
        JET.test_package(Polynomials, ignored_modules=(AnyFrameModule(LinearAlgebra), AnyFrameModule(Base)))
    end
end
