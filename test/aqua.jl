using Aqua
using JET

# The Downgrade workflow promotes the test-only packages into the [deps] of
# Project.toml, which stale_deps and persistent_tasks then trip over. Both are
# still run by the regular CI.
const downgrade_ci = get(ENV, "POLYNOMIALS_DOWNGRADE_CI", "false") == "true"

@testset "Aqua" begin
    Aqua.test_all(Polynomials; stale_deps = !downgrade_ci, persistent_tasks = !downgrade_ci)
end

# Any two-argument `Base.promote` method, whatever its signature, invalidates all
# compiled `promote(a, b)` call sites whose arguments are inferred abstractly.
# Promotion should go through `promote_rule` and `convert` instead.
@testset "No two-argument promote methods" begin
    ms = [m for m in methods(promote)
              if Base.moduleroot(m.module) === Polynomials && m.nargs == 3 && !m.isva]
    @test isempty(ms)
end

if VERSION ≥ v"1.12.0"
    @testset "JET" begin
        JET.test_package(Polynomials, ignored_modules=(AnyFrameModule(LinearAlgebra), AnyFrameModule(Base)))
    end
end
