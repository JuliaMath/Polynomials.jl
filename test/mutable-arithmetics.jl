import MutableArithmetics
const MA = MutableArithmetics

function alloc_test(f, n)
    f() # compile
    @test n == @allocated f()
end


@testset "PolynomialsMutableArithmetics.jl" begin
    d = m = n = 4
    p(d) = Polynomial(big.(1:d))
    z(d) = Polynomial([zero(BigInt) for i in 1:d])
    A = [p(d) for i in 1:m, j in 1:n]
    b = [p(d) for i in 1:n]
    c = [z(2d - 1) for i in 1:m]
    buffer = MA.buffer_for(MA.add_mul, typeof(c), typeof(A), typeof(b))
    @test buffer isa BigInt
    c = [z(2d - 1) for i in 1:m]
    MA.buffered_operate!(buffer, MA.add_mul, c, A, b)
    @test c == A * b
    @test c == MA.operate(*, A, b)
    @test 0 == @allocated MA.buffered_operate!(buffer, MA.add_mul, c, A, b)
end
