using LinearAlgebra
using OffsetArrays, StaticArrays
import Polynomials: indeterminate
import Polynomials: ImmutableDensePolynomial, StandardBasis,MutableSparsePolynomial, MutableDensePolynomial

## Test standard basis polynomials with (nearly) the same tests

#   compare upto trailing  zeros
function  upto_tz(as, bs)
    n,m = findlast.(!iszero, (as,bs))
    n == m || return false
    isnothing(n) &&  return true
    for i in 1:n
        !(as[i] ≈ bs[i]) && return false
    end
    true
end

upto_z(as, bs) = upto_tz(filter(!iszero,as), filter(!iszero,bs))

# compare upto trailing zeros infix operator
==ᵗ⁰(a,b) = upto_tz(a,b)
==ᵗᶻ(a,b) = upto_z(a,b)

==ᵟ(a,b) = (a == b)
==ᵟ(a::FactoredPolynomial, b::FactoredPolynomial) = a ≈ b

_isimmutable(::Type{P}) where {P <: Union{ImmutablePolynomial, FactoredPolynomial, ImmutableDensePolynomial{StandardBasis}}} = true
_isimmutable(P) = false


Ps = (ImmutablePolynomial, Polynomial, SparsePolynomial, LaurentPolynomial, FactoredPolynomial,
      MutableDensePolynomial{StandardBasis},ImmutableDensePolynomial{StandardBasis}, MutableSparsePolynomial{StandardBasis}
      )

@testset "Construction" begin
    @testset for coeff in Any[
          Int64[1, 1, 1, 1],
          Float32[1, -4, 2],
          ComplexF64[1 - 1im, 2 + 3im],
          [3 // 4, -2 // 1, 1 // 1]
         ]

        @testset for P in Ps
            p = P(coeff)
            @test coeffs(p) ==ᵗ⁰ coeff
            @test degree(p) == length(coeff) - 1
            @test indeterminate(p) == :x
            P == Polynomial && @test length(p) == length(coeff)
            P == Polynomial && @test size(p) == size(coeff)
            P == Polynomial && @test size(p, 1) == size(coeff, 1)
            P == Polynomial && @test typeof(p).parameters[1] == eltype(coeff)
            if !(eltype(coeff) <: Real && P == FactoredPolynomial) # roots may be complex
                @test eltype(p) == eltype(coeff)
            end
            @test all([-200, -0.3, 1, 48.2] .∈ Polynomials.domain(p))

            ## issue #316
            @test_throws InexactError P{Int,:x}([1+im, 1])
            @test_throws InexactError P{Int}([1+im, 1], :x)
            @test_throws InexactError P{Int,:x}(1+im)
            @test_throws InexactError P{Int}(1+im)

            ## issue #395
            v = [1,2,3]
            @test P(v) == P(v,:x) == P(v,'x') == P(v,"x") == P(v, Polynomials.Var(:x))

            ## issue #452
            ps = (1, 1.0)
            P != FactoredPolynomial && @test eltype(P(ps)) == eltype(promote(ps...))
            ## issue 464
            @variable z
            @test z^2 + 2 == Polynomial([2,0,1], :z)

            ## issue 457
            @test (@test_deprecated Polynomials.order(p)) == length(coeff) - 1
        end
    end
end

@testset "Mapdomain" begin
    @testset for P in Ps
        x = -30:20
        mx = mapdomain(P, x)
        @test mx == x

        x = 0.5:0.01:0.6
        mx = mapdomain(P, x)
        @test mx == x
    end
end

# Custom offset vector type to test constructors
struct ZVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    x :: A
    offset :: Int
    function ZVector(x::AbstractVector)
        offset = firstindex(x)
        new{eltype(x),typeof(x)}(x, offset)
    end
end
Base.parent(z::ZVector) = z.x
Base.size(z::ZVector) = size(parent(z))
Base.axes(z::ZVector) = (OffsetArrays.IdentityUnitRange(0:size(z,1)-1),)
Base.getindex(z::ZVector, I::Int) = parent(z)[I + z.offset]

@testset "Other Construction" begin
    @testset for P in Ps

        # Leading 0s
        p = P([1, 2, 0, 0])
        @test coeffs(p) ==ᵗ⁰ [1, 2]
        P == Polynomial && @test length(p) == 2

        # different type
        p = P{Float64}(ones(Int32, 4))
        @test coeffs(p) ==ᵗ⁰ ones(Float64, 4)

        p = P(30)
        @test coeffs(p) ==ᵗ⁰ [30]

        p = zero(P{Int})
        @test coeffs(p) ==ᵗ⁰ [0]

        p = one(P{Int})
        @test coeffs(p) ==ᵗ⁰ [1]

        pNULL = P(Int[])
        @test iszero(pNULL)
        P != LaurentPolynomial && @test degree(pNULL) == -1

        p0 = P([0])
        @test iszero(p0)
        P != LaurentPolynomial && @test degree(p0) == -1

        # P(2) is  2 (not  2p₀)  convert(Polynomial, P(s::Number)) = Polynomial(s)
        @test convert(Polynomial, P(2)) ≈ Polynomial(2)
        @test P(2)  ≈ 2*one(P)

        # variable(), P() to generate `x` in given basis
        @test degree(variable(P)) == 1
        @test variable(P)(1) == 1
        @test degree(P()) == 1
        @test P()(1) == 1
        @test variable(P, :y) == P(:y)

        # test degree, isconstant
        P != LaurentPolynomial &&  @test degree(zero(P)) ==  -1
        @test degree(one(P)) == 0
        @test degree(P(1))  ==  0
        @test degree(P([1]))  ==  0
        @test degree(P(:x)) ==  1
        @test degree(variable(P)) == 1
        @test degree(Polynomials.basis(P,5)) == 5
        @test Polynomials.isconstant(P(1))
        @test !Polynomials.isconstant(variable(P))

    end


end

@testset "Non-number type" begin
    conv = Polynomials.conv
    @testset "T=Polynomial{Int,:y}" begin
        @testset for P in (Polynomial,)

            T = P{Int, :y}
            a,b,c = T([1]), T([1,2]), T([1,2,3])
            p = P([a,b,c])
            q = P([a,b])
            s = 2
            d = c

            # scalar product
            @test s*p == P([s*cᵢ for cᵢ ∈ [a,b,c]])
            @test p*s == P([cᵢ*s for cᵢ ∈ [a,b,c]])
            @test_throws ArgumentError d*p == P([d*cᵢ for cᵢ ∈ [a,b,c]]) # can't fix
            @test_throws ArgumentError p*d == P([cᵢ*d for cᵢ ∈ [a,b,c]]) # can't fix

            # poly add
            @test +p == p
            @test p + q == P([a+a,b+b,c])
            @test p - q == P([a-a,b-b,c])
            @test p - p == P([0*a])

            # poly mult
            @test p * q == P(conv([a,b,c], [a,b]))
            @test q * p == P(conv([a,b], [a,b, c]))

            # poly powers
            @test p^2 == p * p

            # evaluation
            @test p(s) == a + b * s + c * s * s
            @test p(c) == a + b * c + c * c * c

            # ∂, ∫
            @test derivative(p) == P([b, 2c])
            @test integrate(p) == P([0*a, a, b/2, c/3])

            # matrix element
            pq = [p q]
            @test pq[1] == p
            @test pq[2] == q

            # implicit promotion
            @test p + s == P([a+s, b, c])
            @test_throws Union{ArgumentError, MethodError} p + d == P([a+d, b, c]) # can't fix
            @test p + P([d]) == P([a+d,b,c])

            ps = [p s]
            @test ps[1] == p
            @test ps[2] == s
        end
    end

    @testset "T=Matrix (2x2)" begin
        @testset for P ∈ (Polynomial, ImmutablePolynomial)
            a,b,c = [1 0; 1 1], [1 0; 2 1], [1 0; 3 1]
            p = P([a,b,c])
            q = P([a,b])
            s = 2
            d = [4 1; 1 0]

            # scalar product
            @test s*p == P([s*cᵢ for cᵢ ∈ [a,b,c]])
            @test p*s == P([cᵢ*s for cᵢ ∈ [a,b,c]])
            @test d*p == P([d*cᵢ for cᵢ ∈ [a,b,c]])
            @test p*d == P([cᵢ*d for cᵢ ∈ [a,b,c]])

            # poly add
            @test +p == p
            @test p + q == P([a+a,b+b,c])
            @test p - q == P([a-a,b-b,c])

            # poly mult
            @test p * q == P(conv([a,b,c], [a,b]))
            @test q * p == P(conv([a,b], [a,b, c]))

            # poly powers
            @test p^2 == p * p

            # evaluation
            @test p(s) == a + b * s + c * s * s
            @test p(c) == a + b * c + c * c * c

            # ∂, ∫
            @test derivative(p) == P([b, 2c])
            @test integrate(p) == P([0*a, a, b/2, c/3])

            # matrix element
            @test [p q][1] == p
            @test [p q][2] == q

            # implicit promotion
            @test_throws MethodError p + s == P([a+s, b, c]) # OK, no a + s
            @test p + d == P([a+d, b, c])
            @test p + P([d]) == P([a+d,b,c])

            @test_throws MethodError [p s][1] == p # no promotion T(s)
            @test_throws MethodError [p s][2] == s
        end
    end


    @testset "T=Vector{Int}" begin
        @testset for P ∈ (Polynomial, ImmutablePolynomial)
            a,b,c = [1,0,0], [1,1,0], [1,1,1]
            p = P([a,b,c])
            q = P([a,b])
            s = 2
            d = [1,2,3]

            # scalar product
            @test s*p == P([s*cᵢ for cᵢ ∈ [a,b,c]])
            @test p*s == P([cᵢ*s for cᵢ ∈ [a,b,c]])
            @test_throws MethodError d*p == P([d*cᵢ for cᵢ ∈ [a,b,c]]) # Ok, no * for T
            @test_throws MethodError p*d == P([cᵢ*d for cᵢ ∈ [a,b,c]]) # Ok, no * for T

            # poly add
            @test +p == p
            @test p + q == P([a+a,b+b,c])
            @test p - q == P([a-a,b-b,c])

            # poly mult
            @test_throws MethodError p * q == P(conv([a,b,c], [a,b])) # Ok, no * for T
            @test_throws MethodError q * p == P(conv([a,b], [a,b, c])) # Ok, no * for T

            # poly powers
            @test_throws MethodError p^2  # Ok, no * for T
            @test_throws MethodError p * p

            # evaluation
            @test p(s) == a + b * s + c * s * s
            @test_throws MethodError p(c) == a + b * c + c * c * c # OK, no b * c

            # ∂, ∫
            @test derivative(p) == P([b, 2c])
            @test integrate(p) == P([0*a, a, b/2, c/3])


            # matrix element
            @test [p q][1] == p
            @test [p q][2] == q

            # implicit promotion
            @test_throws MethodError p + s == P([a+s, b, c])  # OK, no a + s
            @test  p + d == P([a+d, b, c])
            @test p + P([d]) == P([a+d,b,c])

            @test_throws MethodError [p s][1] == p # no promotion T(s)
            @test_throws MethodError [p s][2] == s
        end
    end


    # eval(quote
    #      using StaticArrays
    #      end)
        @testset "T=SA" begin
            @testset for P ∈ (Polynomial, ImmutablePolynomial )
                a,b,c = SA[1 0; 1 1], SA[1 0; 2 1], SA[1 0; 3 1]
                p = P([a,b,c])
                q = P([a,b])
                s = 2
                d = SA[4 1; 1 0]

                # scalar product
                @test s*p == P([s*cᵢ for cᵢ ∈ [a,b,c]])
                @test p*s == P([cᵢ*s for cᵢ ∈ [a,b,c]])
                @test d*p == P([d*cᵢ for cᵢ ∈ [a,b,c]])
                @test p*d == P([cᵢ*d for cᵢ ∈ [a,b,c]])

                # poly add
                @test +p == p
                @test p + q == P([a+a,b+b,c])
                @test p - p == P([0*a])

                # poly mult
                @test p * q == P(conv([a,b,c], [a,b]))
                @test q * p == P(conv([a,b], [a,b, c]))

                # poly powers
                @test p^2 == p * p


                # evaluation
                @test p(s) == a + b * s + c * s * s
                @test p(c) == a + b * c + c * c * c

                # ∂, ∫
                @test derivative(p) == P([b, 2c])
                @test integrate(p) == P([0*a, a, b/2, c/3])

                # matrix element
                @test [p q][1] == p
                @test [p q][2] == q

                # implicit promotion
                # @test_broken p + s == P([a .+ s, b, c]) # should error, doesn't
                @test p + d == P([a + d, b, c])
                @test p + P([d]) == P([a + d,b,c])

                @test_throws MethodError [p s][1] == p # no promotion T(s)
                @test_throws MethodError [p s][2] == s #
            end
        end
end

@testset "OffsetVector" begin
    as = ones(3:4)
    bs = parent(as)


    @testset for P in Ps
        # LaurentPolynomial accepts OffsetArrays; others throw warning
        if P ∈ (LaurentPolynomial, MutableDensePolynomial{StandardBasis})
            @test LaurentPolynomial(as) == LaurentPolynomial(bs, 3)
        else
            @test P(as) == P(bs)
            @test P{eltype(as)}(as) == P{eltype(as)}(bs)
            # (Or throw an error?)
            # @test_throws ArgumentError P(as)
            # @test P{eltype(as)}(as) == P{eltype(as)}(bs)
        end
    end

    a = [1,1]
    b = OffsetVector(a, axes(a))
    c = ZVector(a)
    d = ZVector(b)
    @testset for P in Ps
        if P == LaurentPolynomial && continue
            @test P(a) == P(b) == P(c) == P(d)
        end

    end
end


@testset "Arithmetic" begin

    @testset for P in Ps
        pNULL = P(Int[])
        p0 = P([0])
        p1 = P([1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        p2 = P([1,1,0,0])
        p3 = P([1,2,1,0,0,0,0])
        p4 = P([1,3,3,1,0,0])
        p5 = P([1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        pN = P([276,3,87,15,24,0])
        pR = P([3 // 4, -2 // 1, 1 // 1])

        # type stability of the default constructor without variable name
        if !(P ∈ (LaurentPolynomial, ImmutablePolynomial, FactoredPolynomial, ImmutableDensePolynomial{StandardBasis}))
            @inferred P([1, 2, 3])
            @inferred P([1,2,3], Polynomials.Var(:x))
        end

        @test p3 == P([1,2,1])
        @test pN * 10 == P([2760, 30, 870, 150, 240])
        @test pN / 10.0 ≈ P([27.6, 0.3, 8.7, 1.5, 2.4])
        @test 10 * pNULL + pN ==ᵟ pN
        @test 10 * p0 + pN ==ᵟ pN
        @test p5 + 2 * p1 == P([3,4,6,4,1])
        @test 10 * pNULL - pN ==ᵟ -pN
        @test p0 - pN ==ᵟ -pN
        @test p5 - 2 * p1 == P([-1,4,6,4,1])
        @test p2 * p2 * p2 == p4
        @test p2^4 == p5
        @test pNULL^3 == pNULL
        @test pNULL * pNULL == pNULL

        if P === Polynomial
            # type stability of multiplication
            @inferred 10 * pNULL
            @inferred 10 * p0
            @inferred p2 * p2
            @inferred p2 * p2
        end

        @test pNULL + 2 == p0 + 2 == 2 + p0 == P([2])
        @test p2 - 2 == -2 + p2 == P([-1,1])
        @test 2 - p2 == P([1,-1])

    end

    @testset "generic arithmetics" begin
        P = Polynomial
        # define a set algebra
        struct Setalg  # not a number
            content::Vector{Int}
        end
        Base.:(+)(a::Setalg, b::Setalg) = Setalg(a.content ∪ b.content)
        Base.:(*)(a::Setalg, b::Setalg) = Setalg(vec([x * y for x in a.content, y in b.content]))
        Base.zero(::Setalg) = Setalg(Int[])
        Base.one(::Setalg) = Setalg(Int[1])
        Base.zero(::Type{Setalg}) = Setalg(Int[])
        Base.one(::Type{Setalg}) = Setalg(Int[1])
        Base.:(==)(a::Setalg, b::Setalg) = a.content == b.content

        a = Setalg([1])
        b = Setalg([4,2])
        pNULL = P(Setalg[])
        p0 = P([a])
        p1 = P([a, b, b])
        @test pNULL * p0 == pNULL
        @test pNULL * p1 == pNULL
        @test p0 * p1 == p1
    end

    @testset for P in Ps # ensure promotion of scalar +,*,/
        p = P([1,2,3])
        @test p + 0.5 ==ᵟ P([1.5, 2.0, 3.0])
        @test p / 2  == P([1/2, 1.0, 3/2])
        @test p * 0.5 == P([1/2, 1.0, 3/2])
    end

    # ensure  promotion of +,*; issue 215
    @testset for P in Ps
        p,q = P([1,2,3]), P(im, :θ)
        @test p+q == P([1+im, 2, 3])
        @test p*q ==ᵟ P(im*[1,2,3])
    end

    # Laurent polynomials and scalar operations
    cs = [1,2,3,4]
    p = LaurentPolynomial(cs, -3)
    @test p*3 == LaurentPolynomial(cs .* 3, -3)
    @test 3*p == LaurentPolynomial(3 .* cs, -3)

    # LaurentPolynomial has an inverse for monomials
    x = variable(LaurentPolynomial)
    @test Polynomials.isconstant(x * inv(x))
    @test_throws ArgumentError inv(x + x^2)

    # issue #395
    @testset for P ∈ Ps
        P ∈ (FactoredPolynomial, ImmutablePolynomial, ImmutableDensePolynomial{StandardBasis} ) && continue
        p = P([2,1], :s)
        @inferred -p # issue #395
        @inferred 2p
        @inferred p + p
        @inferred p * p
        @inferred p^3
    end

    # evaluation at special cases (degree 0,1; evaluate at 0)
    @testset for P ∈ Ps
        for T ∈ (Int, Float16, Float64, Complex{Float64})
            p₀ = zero(P{T,:X})
            @test p₀(0) == zero(T) == Polynomials.constantterm(p₀)
            @test p₀(1) == zero(T)
            p₁ = P(T[2])
            @test p₁(0) == T(2) == Polynomials.constantterm(p₁)
            @test p₁(1) == T(2)
        end
    end

    # evaluation at special cases different number types
    @testset for P ∈ Ps
        P ∈ (SparsePolynomial, FactoredPolynomial, MutableSparsePolynomial{StandardBasis}) && continue
        # vector coefficients
        v₀, v₁ = [1,1,1], [1,2,3]
        p₁ = P([v₀])
        @test p₁(0) == v₀  == Polynomials.constantterm(p₁)
        @test_throws MethodError (0 * p₁)(0) # no zero(Vector{Int})
        p₂ = P([v₀, v₁])
        @test p₂(0) == v₀ == Polynomials.constantterm(p₂)
        @test p₂(2) == v₀ + 2v₁

        # matrix arguments
        # for matrices like pₒI + p₁X + p₂X² + ⋯
        p = P([1])
        x = [1 2; 3 4]
        @test p(x) == 1*I
        @test (0p)(x) == 0*I
        p = P([1,2])
        @test p(x) == 1*I + 2*x
        p = P([1,2,3])
        @test p(x) == 1*I + 2*x + 3x^2
    end

    # p - p requires a zero
    @testset for P ∈ Ps
        P ∈ (LaurentPolynomial, SparsePolynomial,
             FactoredPolynomial, MutableSparsePolynomial{StandardBasis}) && continue
        for v ∈ ([1,2,3],
                 [[1,2,3],[1,2,3]],
                 [[1 2;3 4], [3 4; 5 6]]
                 )
            p = P(v)
            @test p - p == 0*p
        end
    end

    # issue #495, (scalar div fix)
    𝐐 = Rational{Int}
    v = Polynomial{𝐐}([0//1])
    @test eltype(integrate(v)) == 𝐐
end

@testset "Divrem" begin
    @testset for P in  Ps
        P == FactoredPolynomial && continue
        p0 = P([0])
        p1 = P([1])
        p2 = P([5, 6, -3, 2 ,4])
        p3 = P([7, -3, 2, 6])
        p4 = p2 * p3
        pN = P([276,3,87,15,24,0])
        pR = P([3 // 4, -2 // 1, 1 // 1])

        @test all(divrem(p4, p2) .==ᵟ (p3, zero(p3)))
        @test p3 % p2 ==ᵟ p3
        @test all((map(abs, coeffs(p2 ÷ p3 - P([1 / 9,2 / 3])))) .< eps())
        @test all(divrem(p0, p1) .==ᵟ (p0, p0))
        @test all(divrem(p1, p1) .==ᵟ (p1, p0))
        @test all(divrem(p2, p2) .==ᵟ (p1, p0))
        @test all(divrem(pR, pR) .==ᵟ (one(pR), zero(pR)))
        @test_throws DivideError p1 ÷ p0
        @test_throws DivideError divrem(p0, p0)

        # issue #235
        num = P([0.8581454436924945, 0.249671302254737, 0.8048498901050951, 0.1922713965697087]) # degree 3 polynomial
        den = P([0.9261520696359462, 0.07141031902098072, 0.378071465860349]) # degree 2 polynomial
        q, r = divrem(num,den)  # expected degrees: degree(q) = degree(num)-degree(den) = 1, degree(r) = degree(den)-1 = 1
        @test num ≈ den*q+r  # true
        @test degree(q) == 1 # true
        degree(r) < degree(den)
    end
end

@testset "Comparisons" begin
    @testset for P in Ps
        pX = P([1, 2, 3, 4, 5])
        pS1 = P([1, 2, 3, 4, 5], "s")
        pS2 = P([1, 2, 3, 4, 5], 's')
        pS3 = P([1, 2, 3, 4, 5], :s)
        @test pX != pS1
        @test pS1 == pS2
        @test pS1 == pS3

        @test indeterminate(pS1 + pS1) == indeterminate(pS1)
        @test indeterminate(pS1 - pS1) == indeterminate(pS1)
        @test indeterminate(pS1 * pS1) == indeterminate(pS1)
        @test indeterminate(pS1 ÷ pS1) == indeterminate(pS1)
        @test indeterminate(pS1 % pS1) == indeterminate(pS1)

        @test_throws ArgumentError pS1 + pX
        @test_throws ArgumentError pS1 - pX
        @test_throws ArgumentError pS1 * pX
        @test_throws ArgumentError pS1 ÷ pX
        @test_throws ArgumentError pS1 % pX

        # Testing copying.
        pcpy1 = P([1,2,3,4,5], :y)
        pcpy2 = copy(pcpy1)
        @test pcpy1 == pcpy2

        # Check for isequal
        p1 = P([1.0, -0.0, 5.0, Inf])
        p2 = P([1.0,  0.0, 5.0, Inf])
        !(P ∈ (FactoredPolynomial, SparsePolynomial, MutableSparsePolynomial{StandardBasis})) && (@test p1 == p2 && !isequal(p1, p2))  # SparsePolynomial doesn't store -0.0,  0.0.

        p3 = P([0, NaN])
        @test p3 === p3 && p3 ≠ p3 && isequal(p3, p3)

        p = fromroots(P, [1,2,3])
        q = fromroots(P, [1,2,3])
        @test hash(p) == hash(q)

        p1s = P([1,2], :s)
        p1x = P([1,2], :x)
        p2s = P([1], :s)

        @test p1s == p1s
        @test p1s ≠ p1x
        @test p1s ≠ p2s

        @test_throws ArgumentError p1s ≈ p1x
        @test p1s ≉ p2s
        @test p1s ≈ P([1,2.], :s)

        @test p2s ≈ 1.0 ≈ p2s
        @test p2s == 1.0 == p2s
        @test p2s ≠ 2.0 ≠ p2s
        @test p1s ≠ 2.0 ≠ p1s

        @test nnz(map(P, sparse(1.0I, 5, 5))) == 5

        @test P([0.5]) + 2 == P([2.5])
        @test 2 - P([0.5]) == P([1.5])

        # check ≈ for P matches usage for Vector{T} (possibly padded with trailing zeros)
        @test (P([NaN]) ≈ P([NaN]))               == ([NaN] ≈ [NaN]) # false
        @test (P([NaN]) ≈ NaN)                    == (false)
        @test (P([Inf]) ≈ P([Inf]))               == ([Inf] ≈ [Inf]) # true
        @test (P([Inf]) ≈ Inf)                    == (true)
        !(P <: FactoredPolynomial) && @test (P([1,Inf]) ≈ P([0,Inf])) == ([1,Inf] ≈ [0,Inf]) # false
        !(P <: FactoredPolynomial) && @test (P([1,NaN,Inf]) ≈ P([0,NaN, Inf])) == ([1,NaN,Inf] ≈ [0,NaN, Inf]) #false
        @test (P([eps(), eps()]) ≈ P([0,0]))      == ([eps(), eps()] ≈ [0,0]) # false
        @test (P([1,eps(), 1]) ≈ P([1,0,1]))      == ([1,eps(), 1] ≈ [1,0,1]) # true
        !(P <: FactoredPolynomial) && @test (P([1,2]) ≈ P([1,2,eps()])) == ([1,2,0] ≈ [1,2,eps()])

        # NaN poisons comparison
        @test !(P([NaN, 1.0, 2.0]) ≈ P([NaN, 1.0, 2.0]))

        # check how ==, ===, isapprox ignore variable mismatch when constants are involved, issue #217, issue #219
        @test zero(P, :x) == zero(P, :y)
        @test one(P, :x) == one(P, :y)
        @test !(variable(P, :x) == variable(P,:y))

        @test !(zero(P, :x) === zero(P, :y))
        @test !(one(P, :x) === one(P, :y))
        @test !(variable(P, :x) === variable(P,:y))

        @test zero(P, :x) ≈ zero(P, :y)
        @test one(P, :x) ≈ one(P, :y)
        @test (variable(P, :x) ≈ variable(P, :x))
        @test_throws ArgumentError variable(P, :x) ≈ variable(P, :y)

    end

    ## Issue 408
    p = Polynomial([1,2,3])
    q = ChebyshevT([1,2,3])
    @test p != q
end

@testset "Fitting" begin
    @testset for P in Ps
        P <: FactoredPolynomial && continue
        xs = range(0, stop = π, length = 10)
        ys = sin.(xs)

        p = fit(P, xs, ys)
        y_fit = p.(xs)
        abs_error = abs.(y_fit .- ys)
        @test maximum(abs_error) <= 0.03

        p = fit(P, xs, ys, 2)
        y_fit = p.(xs)
        abs_error = abs.(y_fit .- ys)
        @test maximum(abs_error) <= 0.03

        # Test weighted
        @testset for W in [1, ones(size(xs)), diagm(0 => ones(size(xs)))]
            p = fit(P, xs, ys, 2, weights = W)
            @test p.(xs) ≈ y_fit
        end


        # Getting error on passing Real arrays to polyfit #146
        xx = Real[20.0, 30.0, 40.0]
        yy = Real[15.7696, 21.4851, 28.2463]
        fit(P, xx, yy, 2)

        # issue #214 --  should error
        @test_throws ArgumentError fit(Polynomial, rand(2,2), rand(2,2))

        # issue #268 -- inexacterror
        @test fit(P, 1:4, 1:4, var=:x) ≈ variable(P{Float64}, :x)
        @test fit(P, 1:4, 1:4, 1, var=:x) ≈ variable(P{Float64}, :x)

        # issue #467, fit specific degrees only
        p = fit(P, xs, ys, 1:2:9)
        @test norm(p.(xs) - ys) ≤ 1e-4

        # issue 467: with constants
        p = fit(P, xs, ys, 3:2:9, Dict(1 => 1))
        @test norm(p.(xs) - ys) ≤ 1e-3

    end



    f(x) = 1/(1 + 25x^2)
    N = 250; xs = [cos(j*pi/N) for j in N:-1:0];
    q = fit(ArnoldiFit, xs, f.(xs));
    @test maximum(abs, q(x) - f(x) for x ∈ range(-1,stop=1,length=500)) < 10eps()
    q = fit(ArnoldiFit, xs, f.(xs), 100);
    @test maximum(abs, q(x) - f(x) for x ∈ range(-1,stop=1,length=500)) < sqrt(eps())


    # test default   (issue  #228)
    fit(1:3,  rand(3))

    # weight test (PR #291)
    # we specify w^2.
    x = range(0, stop=pi, length=30)
    y = sin.(x)
    wts = 1 ./ sqrt.(1 .+ x)
    # cs = numpy.polynomial.polynomial.polyfit(x, y, 4, w=wts)
    cs = [0.0006441172319036863, 0.985961582190304, 0.04999233434370933, -0.23162369757680354, 0.036864056406570644]
    @test maximum(abs, coeffs(fit(x, y, 4, weights=wts.^2)) - cs) <= sqrt(eps())
end

@testset "Values" begin
    @testset for P in Ps
        pNULL = P(Int[])
        p0 = P([0])
        p1 = P([1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        p2 = P([1,1,0,0])
        p3 = P([1,2,1,0,0,0,0])
        p4 = P([1,3,3,1,0,0])
        p5 = P([1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        pN = P([276,3,87,15,24,0])
        pR = P([3 // 4, -2 // 1, 1 // 1])


        @test fromroots(P, Int[])(2.) == 1.
        @test pN(-.125) ≈ 276.9609375
        @test pNULL(10) == 0
        @test p0(-10) == 0
        @test fromroots(P, [1 // 2, 3 // 2])(1 // 2) == 0 // 1

        # Check for Inf/NaN operations
        p1 = P([Inf, Inf])
        p2 = P([0, Inf])
        @test p1(Inf) == Inf
        if !(P <: FactoredPolynomial)
            @test isnan(p1(-Inf))
            @test isnan(p1(0))
            @test p2(-Inf) == -Inf
        end

        # issue #189
        p = P([0,1,2,3])
        A = [0 1; 0  0];
        @test  p(A) ==ᵟ A + 2A^2 + 3A^3

        # issue #209
        ps  = [P([0,1]), P([0,0,1])]
        @test Polynomials.evalpoly.(1/2, ps) ≈ [p(1/2)  for  p  in ps]

    end


    # constant polynomials and type
    Ts = (Int, Float32, Float64, Complex{Int}, Complex{Float64})
    @testset for P in (Polynomial, ImmutablePolynomial, SparsePolynomial)
        @testset for T in Ts
            @testset for S in Ts
                c = 2
                p = P{T}(c)
                x = one(S)
                y = p(x)
                @test y === c * one(T)*one(S)
                q = P{T}([c,c])
                @test typeof(q(x)) == typeof(p(x))
            end
        end
    end

    @testset for P in Ps
        p = P(1)
        x = [1 0; 0 1]
        y = p(x)
        @test y ≈ x

        # Issue #208 and  type of output
        p1=P([1//1])
        p2=P([0, 0.9])
        p3=p1(p2)
        @test isa(p3, P)
        @test eltype(p3) == eltype(p2)
    end

    # compensated_horner
    # polynomial evaluation for polynomials with large condition numbers
    @testset for P in (Polynomial, ImmutablePolynomial, SparsePolynomial)
        x = variable(P{Float64})
        f(x) = (x - 1)^20
        p = f(x)
        e₁ = abs( (f(4/3) - p(4/3))/ p(4/3) )
        e₂ = abs( (f(4/3) - Polynomials.compensated_horner(coeffs(p), 4/3))/ p(4/3) )
        λ = cond(p, 4/3)
        u = eps()/2
        @test λ > 1/u
        @test e₁ <= 2 * 20 * u * λ
        @test e₁ > u^(1/4)
        @test e₂ <= u + u^2 * λ * 100
    end
end

@testset "Conversion" begin

    X = :x
    @testset for P in Ps
        if !(P ∈ (ImmutablePolynomial,  ImmutableDensePolynomial{StandardBasis}))
            p = P([0,one(Float64)])
            @test P{Complex{Float64},X} == typeof(p + 1im)
            @test P{Complex{Float64},X} == typeof(1im - p)
            @test P{Complex{Float64},X} == typeof(p * 1im)
        else
            p = P([0,one(Float64)])
            N=2
            @test P{Complex{Float64},X,N} == typeof(p + 1im)
            @test P{Complex{Float64},X,N} == typeof(1im - p)
            @test P{Complex{Float64},X,N} == typeof(p * 1im)
        end
    end

    # unnecessary copy in convert #65
    p1 = Polynomial([1,2])
    p2 = convert(Polynomial{Int}, p1)
    p2[3] = 3
    @test p1[3] == 3

    # issue #287
    p = LaurentPolynomial([1], -5)
    @test p ≈ convert(LaurentPolynomial{Float64}, p)

    # issue #358 `P(p::AbstractPolynomial)` should be `convert(P, p)` not `P(pᵢ for pᵢ ∈ p))`
    x² = Polynomial([0,0,1], :x)
    @testset for P ∈ (ImmutablePolynomial, SparsePolynomial, ChebyshevT)
        @test P(x²) == convert(P, x²)
        Q = P{Float64}
        @test Q(x²) == convert(Q, x²)
    end

    # preserve eltype in SparsePolynomial
    s = SparsePolynomial(Dict(1=>3, 2=>4))
    s2 = SparsePolynomial(s)
    @test s2 isa typeof(s)
    @test s2 == s
    s3 = SparsePolynomial{Float64}(s)
    @test s3 isa SparsePolynomial{Float64,indeterminate(s)}

    # conversions between pairs of polynomial types
    c = [1:5;]
    Psexact = (ImmutablePolynomial, Polynomial, SparsePolynomial, LaurentPolynomial)
    @testset for P1 in Ps
        p = P1(c)
        @testset for P2 in Psexact
            @test convert(P2, p) == p
        end
        @test convert(FactoredPolynomial, p) ≈ p
    end

    # reinterpret coefficients
    for P in (ImmutablePolynomial, Polynomial, SparsePolynomial, LaurentPolynomial, FactoredPolynomial)
        for T in (Float64, Rational)
            xs = [1,2,3]
            p = fromroots(P,xs)
            @test Polynomials.copy_with_eltype(T, p) == fromroots(P, T.(xs))
            @test Polynomials.copy_with_eltype(T, Val(:u), p) == fromroots(P, T.(xs); var=:u)
            P == ImmutablePolynomial && continue
            @inferred Polynomials.copy_with_eltype(T, Val(:u), p)
            @inferred Polynomials.copy_with_eltype(T, p)
        end
    end
end

@testset "Roots" begin
    @testset for P in Ps

        pNULL = P(Int[])
        p0 = P([0])
        p1 = P([1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        p2 = P([1,1,0,0])
        p3 = P([1,2,1,0,0,0,0])
        p4 = P([1,3,3,1,0,0])
        p5 = P([1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        pN = P([276,3,87,15,24,0])
        pR = P([3 // 4, -2 // 1, 1 // 1])

        # From roots
        r = [2, 3]
        @test fromroots(r) ==ᵟ Polynomial([6, -5, 1])
        p = fromroots(P, r)
        @test p ==ᵟ P([6, -5, 1])
        @test sort(roots(p)) ≈ r

        @test roots(p0) == roots(p1) == roots(pNULL) == []
        @test eltype(roots(p0)) == eltype(roots(p1)) == eltype(roots(pNULL))
        !(P <: FactoredPolynomial) && @test eltype(roots(pNULL)) == Float64
        @test P == LaurentPolynomial ? roots(variable(P)) == [0.0] : roots(P([0,1,0])) == [0.0]

        @test roots(p2) == [-1]
        if !(P <: FactoredPolynomial)
            a_roots = [c for c in coeffs(copy(pN))]
            @test all(map(abs, sort(roots(fromroots(a_roots))) - sort(a_roots)) .< 1e6)
        end
        @test length(roots(p5)) == 4
        @test roots(pNULL) == []
        @test sort(roots(pR)) == [1 // 2, 3 // 2]

        @test sort(roots(LaurentPolynomial([24,10,-15,0,1],-2,:z)))≈[-4.0,-1.0,2.0,3.0]

        A = [1 0; 0 1]
        @test fromroots(A) == Polynomial(Float64[1, -2, 1])
        p = fromroots(P, A)
        @test p == P(Float64[1, -2, 1])
        @test roots(p) ≈ sort!(eigvals(A), rev = true)

        x = variable()
        plarge = 8.362779449448982e41 - 2.510840694154672e57x + 4.2817430781178795e44x^2 - 1.6225927682921337e31x^3 + 1.0x^4  # #120
        @test length(roots(plarge)) == 4

        @test begin
            a = P([1,1,1])*P([1,0.5,1])*P([1,1])    # two complex conjugate pole pairs and one real pole
            r = roots(a)
            b = fromroots(r)
            (b ≈ a) & isreal(coeffs(b))    # the coeff should be real
        end
    end
end

@testset "multroot" begin
    @testset for P in (Polynomial, ImmutablePolynomial, SparsePolynomial)
        rts = [1.0, sqrt(2), sqrt(3)]
        ls = [2, 3, 4]
        x = variable(P{Float64})
        p = prod((x-z)^l for (z,l) in zip(rts, ls))
        out = Polynomials.Multroot.multroot(p)
        @test all(out.values .≈ rts)
        @test all(out.multiplicities .≈ ls)
        @test out.ϵ <= sqrt(eps())
        @test out.κ * out.ϵ < sqrt(eps())  # small forward error
        # one for which the multiplicities are not correctly identified
        n = 4
        q = p^n
        out = Polynomials.Multroot.multroot(q)
        @test (out.multiplicities == n*ls) || (out.κ * out.ϵ > sqrt(eps()))  # large  forward error, l misidentified
        # with right manifold it does yield a small forward error
        zs′ = Polynomials.Multroot.pejorative_root(q, rts .+ 1e-4*rand(3), n*ls)
        @test prod(Polynomials.Multroot.stats(q, zs′, n*ls))  < sqrt(eps())
        # bug with monomial
        T = Float64
        x = variable(P{T})
        out = Polynomials.Multroot.multroot(x^3)
        @test out.values == zeros(T,1)
        @test out.multiplicities == [3]
        # bug with constant
        out = Polynomials.Multroot.multroot(P(1))
        @test isempty(out.values)
        @test isempty(out.multiplicities)
    end
end

@testset "critical points" begin
    for P in (Polynomial, ImmutablePolynomial)
        p = fromroots(P, [-1,-1, 2]) |> integrate
        cps = Polynomials.critical_points(p, (-5,5); endpoints=false)
        @test all(cps .≈ [-1, 2])
        cps = Polynomials.critical_points(p, (0,5); endpoints=false)
        @test all(cps .≈ [2])

        cps = Polynomials.critical_points(p)
        m, i = findmin(p, cps)
        @test m ≈ -6.0
        x = argmin(p, cps)
        @test x ≈ 2.0
        mn, mx = extrema(p, cps)
        @test mn ≈ -6.0 && isinf(mx)
    end
end

@testset "Integrals and Derivatives" begin
    # Integrals derivatives
    @testset for P in Ps
        P <: FactoredPolynomial && continue
        pNULL = P(Int[])
        p0 = P([0])
        p1 = P([1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        p2 = P([1,1,0,0])
        p3 = P([1,2,1,0,0,0,0])
        p4 = P([1,3,3,1,0,0])
        p5 = P([1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
        pN = P([276,3,87,15,24,0])
        pR = P([3 // 4, -2 // 1, 1 // 1])

        c = [1, 2, 3, 4]
        p = P(c)

        der = derivative(p)
        @test coeffs(der) ==ᵗ⁰ [2, 6, 12]
        int = integrate(der, 1)
        @test coeffs(int) ==ᵗ⁰ c


        @test derivative(pR) == P([-2 // 1,2 // 1])
        @test derivative(p3) == P([2,2])
        @test derivative(p1) == derivative(p0) == derivative(pNULL) == pNULL
        @test_throws ArgumentError derivative(pR, -1)
        @test integrate(P([1,1,0,0]), 0, 2) == 4.0

        @test derivative(integrate(pN)) ==ᵟ pN
        @test integrate(pNULL, 1) ==ᵟ p1
        rc = Rational{Int64}[1,2,3]
        @test integrate(P(rc)) == P{eltype(rc)}([0, 1, 1, 1])


        P <: FactoredPolynomial && continue

        @testset for i in 1:10
            p = P(rand(1:5, 6))
            @test degree(truncate(p - integrate(derivative(p)), atol=1e-13)) <= 0
            @test degree(truncate(p - derivative(integrate(p)), atol=1e-13)) <= 0
        end

        # Handling of `NaN`s
        p     = P([NaN, 1, 5])
        pder  = derivative(p)
        pint  = integrate(p)

        @test isnan(p(1)) # p(1) evaluates to NaN
        @test isequal(pder, P([NaN]))
        @test isequal(pint, P([NaN]))

        c = 0.0im
        pint  = integrate(p, c)
        @test isequal(pint, P{promote_type(eltype(p), typeof(c)), :x}([NaN]))

        # Issue with overflow and polyder Issue #159
        @test derivative(P(BigInt[0, 1])^100, 100) == P(factorial(big(100)))
    end

    # Sparse and Laurent test issue #409
    pₛ, pₗ = SparsePolynomial(Dict(-1=>1, 1=>1)), LaurentPolynomial([1,0,1], -1)
    @test pₛ == pₗ && derivative(pₛ) == derivative(pₗ)
    @test_throws ArgumentError integrate(pₛ)
    @test_throws ArgumentError integrate(pₗ)
    qₛ, qₗ = SparsePolynomial(Dict(-2=>1, 1=>1)), LaurentPolynomial([1,0,0,1], -2)
    @test qₛ == qₗ && integrate(qₛ) == integrate(qₗ)
end

@testset "Elementwise Operations" begin
    @testset for P in Ps
        p1  = P([1, 2])
        p2  = P([3, 1.])
        p   = [p1, p2]
        q   = [3, p1]
        if !_isimmutable(P)
            @test q isa Vector{typeof(p1)}
            @test p isa Vector{typeof(p2)}
        else
            @test q isa Vector{<:P{eltype(p1),:x}} # ImmutablePolynomial{Int64,N} where {N}, different  Ns
            @test p isa Vector{<:P{eltype(p2),:x}} # ImmutablePolynomial{Int64,N} where {N}, different  Ns
        end



        psum  = p .+ 3
        pprod = p .* 3
        pmin  = p .- 3
        @test psum  isa Vector{typeof(p2)}
        @test pprod isa Vector{typeof(p2)}
        @test pmin  isa Vector{typeof(p2)}
    end
end

@testset "Chop and Truncate" begin
    # chop and truncate
    @testset for P in Ps
        P <: FactoredPolynomial && continue
        if P == Polynomial
            p = P([1, 1, 1, 1])
            coeffs(p)[end] = 0
            @assert coeffs(p) == [1, 1, 1, 0]
            p = chop(p)
        else
            p = P([1, 1, 1, 0])
        end

        @test coeffs(p) ==ᵗ⁰ [1, 1, 1]
        ## truncation
        p1 = P([1,1] / 10)
        p2 = P([1,2] / 10)
        p3 = P([1,3] / 10)
        psum = p1 + p2 - p3
        @test degree(psum) == 1         # will have wrong degree
        @test degree(truncate(psum)) == 0 # the degree should be correct after truncation

        @test truncate(P([2,1]), rtol = 1 / 2, atol = 0) == P([2])
        @test truncate(P([2,1]), rtol = 1, atol = 0)   == P([0])
        @test truncate(P([2,1]), rtol = 0, atol = 1)   == P([2])

        pchop = P([1, 2, 3, 0, 0, 0])
        pchopped = chop(pchop)
        @test roots(pchop) == roots(pchopped)

    end
end


@testset "As matrix elements" begin

    @testset for P in Ps
        p = P([1,2,3], :x)
        A = [1 p; p^2 p^3]
        @test !issymmetric(A)
        U = A * A'
        @test U[1,2] ≈ U[2,1] # issymmetric with some allowed error for FactoredPolynomial
        P != Polynomials.ImmutableDensePolynomial{Polynomials.StandardBasis} && diagm(0 => [1, p^3], 1=>[p^2], -1=>[p])
    end

    # issue 206 with mixed variable types and promotion
    @testset for P in Ps
        λ = P([0,1],:λ)
        A = [1 λ; λ^2 λ^3]
        P != Polynomials.ImmutableDensePolynomial{Polynomials.StandardBasis} && @test A ==  diagm(0 => [1, λ^3], 1=>[λ], -1=>[λ^2]) # XXX diagm + ImmutableDensePolynomial{StandardBasis} isn't working
        @test all([1 -λ]*[λ^2 λ; λ 1] .== 0)
        @test [λ 1] + [1 λ] == (λ+1) .* [1 1] # (λ+1) not a number, so we broadcast
    end

    # issue 312; using mixed polynomial types within arrays and promotion
    P′ = Polynomial
    r,s = P′([1,2], :x), P′([1,2],:y)
    function _test(x, T,X)
        U = eltype(x)
        Polynomials.constructorof(U) == T && Polynomials.indeterminate(U) == X
    end
    meths = (Base.vect, Base.vcat, Base.hcat)
    @testset for P in (Polynomial, ImmutablePolynomial, SparsePolynomial, LaurentPolynomial)

        p,q = P([1,2], :x), P([1,2], :y)
        P′′ = P == LaurentPolynomial ? P : P′ # different promotion rule

        # * should promote to Polynomial type if mixed (save Laurent Polynomial)
        @testset "promote mixed polys" begin
            @testset for m ∈ meths
                @test _test(m(p,p), P, :x)
                @test _test(m(p,r), P′′, :x)
            end

            @test _test(Base.hvcat((2,1), p, r,[p r]), P′′, :x)

        end

        # * numeric constants should promote to a polynomial, when mixed
        @testset "promote numbers to poly" begin
            @testset for m ∈ meths
                @test _test(m(p,1), P, :x)
                @test _test(m(1,p), P, :x)
                @test _test(m(1,1,p), P, :x)
                @test _test(m(p,1,1), P, :x)
            end

            @test _test(Base.hvcat((3,1), 1, p, r,[1 p r]), P′′, :x)
        end

        # * non-constant polynomials must share the same indeterminate
        @testset "non constant polys share same X" begin
            @testset for m ∈ meths
                @test_throws ArgumentError m(p,q)
                @test_throws ArgumentError m(p,s)
            end

            @test_throws ArgumentError Base.hvcat((2,1), p, q,[p q])
        end


        # * constant polynomials are treated as `P{T,X}`, not elements of `T`
         @testset "constant polys" begin
            @testset for m ∈ meths
                @test _test(m(one(p),1), P, :x)
                @test _test(m(1,one(p)), P, :x)
                @test _test(m(1,1,one(p)), P, :x)
                @test _test(m(one(p),1,1), P, :x)
            end

            @test _test(Base.hvcat((3,1), 1, p, r,[1 p r]), P′′, :x)
        end

        # * Promotion can be forced to mix constant-polynomials
        @testset "Use typed constructor to mix constant polynomials" begin
            𝑷,𝑸 = P{Int,:x}, P{Int,:y} # not typeof(p),... as Immutable carries N
            @test_throws ArgumentError [one(p), one(q)]
            @test eltype(𝑷[one(p), one(q)]) == 𝑷
            @test eltype(𝑸[one(p), one(q)]) == 𝑸
            @test eltype(𝑷[one(p); one(q)]) == 𝑷
            @test eltype(𝑸[one(p); one(q)]) == 𝑸
            @test eltype(𝑷[one(p)  one(q)]) == 𝑷
            @test eltype(𝑸[one(p)  one(q)]) == 𝑸

            @test_throws ArgumentError [1 one(p);
                                        one(p) one(q)]
            @test eltype(𝑷[1 one(p); one(p) one(q)]) == 𝑷
            @test eltype(𝑸[1 one(p); one(p) one(q)]) == 𝑸
        end

        @testset "hvcat" begin
            p,q = P([1,2],:x), P([1,2],:y)

            q1 = [q 1]
            q11 = [one(q) 1]

            @test_throws ArgumentError hvcat((2,1), 1, p, q1)

            @test_throws ArgumentError [1 p; q11]
            @test_throws ArgumentError hvcat((2,1), 1, p, q11)

            𝑷 = P{Int,:x}
            @test eltype(𝑷[1 p; q11]) == 𝑷
            @test eltype(Base.typed_hvcat(𝑷, (2, 1), 1, p, q11)) == 𝑷
        end

    end


end

@testset "Linear Algebra" begin
    @testset for P in Ps
        p = P([3, 4])
        @test norm(p) == 5
        p = P([-1, 3, 5, -2])
        @test norm(p) ≈ 6.244997998398398
        p = P([1 - 1im, 2 - 3im])
        p2 = conj(p)
        @test coeffs(p2) ==ᵗ⁰ [1 + 1im, 2 + 3im]
        @test transpose(p) == p
        !_isimmutable(P) &&  @test transpose!(p) == p
        @test adjoint(Polynomial(im)) == Polynomial(-im) # issue 215
        @test conj(Polynomial(im)) == Polynomial(-im) # issue 215

        @test norm(P([1., 2.])) == norm([1., 2.])
        @test norm(P([1., 2.]), 1) == norm([1., 2.], 1)
    end


    ## Issue #225 and different meanings for "conjugate"
    P = LaurentPolynomial
    p = P(rand(Complex{Float64}, 4), -1)
    z = rand(Complex{Float64})
    s = imag(z)*im
    @test conj(p)(z) ≈ (conj ∘ p ∘ conj)(z)
    @test Polynomials.paraconj(p)(z) ≈ (conj ∘ p ∘ conj ∘ inv)(z)
    @test Polynomials.cconj(p)(s) ≈ (conj ∘ p)(s)

end

@testset "Indexing" begin
    # Indexing
    @testset for P in Ps

        # getindex
        p = P([-1, 3, 5, -2])
        @test p[0] ≈ -1
        @test p[[1, 2]] ==ᵗ⁰ [3, 5]
        @test p[1:2] ==ᵗ⁰ [3, 5]
        @test p[:] ==ᵗ⁰ [-1, 3, 5, -2]

        # setindex
        p1  = P([1,2,1])
        if !_isimmutable(P)
            p1[5] = 1
            @test p1[5] == 1
            @test p1 == P([1,2,1,0,0,1])

            @test p[end] == coeffs(p)[end]

            if P != SparsePolynomial
                p[1] = 2
                @test coeffs(p)[2] == 2
                p[2:3] = [1, 2]
                @test coeffs(p)[3:4] == [1, 2]
                p[0:1] = 0
                @test coeffs(p)[1:2] ==ᵗ⁰ [0, 0]

                p[:] = 1
                @test coeffs(p) ==ᵗ⁰ ones(4)
            end

            p[:] = 0
            @test chop(p) ≈ zero(p)
        end

        p1 = P([1,2,0,3])
        @testset for term in p1
            @test isa(term, eltype(p1))
        end

        !(P <: FactoredPolynomial) && @test eltype(p1) == Int
        @testset for P in Ps
            p1 = P([1,2,0,3])
            @test eltype(collect(p1)) <: eltype(p1)
            @test eltype(collect(Float64, p1)) <: Float64
            @test_throws InexactError collect(Int, P([1.2]))
        end

        p1 = P([1,2,0,3])
        @test length(collect(p1)) == degree(p1) + 1

        @test [p1[idx] for idx in eachindex(p1)] ==ᵗᶻ [1,2,0,3]
    end
end

@testset "Iteration" begin
    p, ip, lp, sp = ps = (Polynomial([1,2,0,4]), ImmutablePolynomial((1,2,0,4)),
                          LaurentPolynomial([1,2,0,4], -2), SparsePolynomial(Dict(0=>1, 1=>2, 3=>4)))
    @testset for pp ∈ ps
        # iteration
        @test all(collect(pp) .== coeffs(pp))

        # keys, values, pairs
        ks, vs, kvs = keys(pp), values(pp), pairs(pp)
        if !isa(pp, SparsePolynomial)
            @test first(ks) == firstindex(pp)
            @test first(vs) == pp[first(ks)]
            @test length(vs) == length(coeffs(pp))
            @test first(kvs) == (first(ks) => first(vs))
        else
            @test first(sort(collect(ks))) == firstindex(pp)
            @test length(vs) == length(pp.coeffs)
        end

        ## monomials
        if !isa(pp, SparsePolynomial)
            i = firstindex(pp)
            @test first(Polynomials.monomials(pp)) == pp[i] * Polynomials.basis(pp,i)
        else
            @test first(Polynomials.monomials(pp)) ∈ [pp[i] * Polynomials.basis(pp,i) for i ∈ keys(pp)]
        end
    end

    @testset "PolynomialKeys/PolynomialValues" begin
        p = Polynomial(Float64[1,2,3])
        pv = Polynomials.PolynomialValues{typeof(p)}(p)
        @test length(pv) == length(p)
        @test size(pv) == size(p)
        @test eltype(pv) == eltype(p)
        pk = Polynomials.PolynomialKeys{typeof(p)}(p)
        @test eltype(pk) == Int
        @test length(pv) == length(p)
        @test size(pv) == size(p)
    end

end


@testset "Copying" begin
    @testset for P in Ps
        pcpy1 = P([1,2,3,4,5], :y)
        pcpy2 = copy(pcpy1)
        @test pcpy1 == pcpy2
    end
end

@testset "GCD" begin
    @testset for P in Ps
        p1 = P([2.,5.,1.])
        p2 = P([1.,2.,3.])

        @test degree(gcd(p1, p2)) == 0          # no common roots
        @test degree(gcd(p1, P(5))) == 0          # ditto
        @test degree(gcd(p1, P(eps(0.)))) == 0          # ditto
        @test degree(gcd(p1, P(0))) == degree(p1) # P(0) has the roots of p1
        @test degree(gcd(p1 + p2 * 170.10734737144486, p2)) == 0          # see, c.f., #122



        p1 = fromroots(P, [1.,2.,3.])
        p2 = fromroots(P, [1.,2.,6.])
        res = roots(gcd(p1, p2))
        @test 1. ∈ res
        @test 2. ∈ res
    end


    # issue 240
    P = Polynomial

    a = P([0.8457170323029561, 0.47175077674705257,  0.9775441940117577]);
    b = P([0.5410010714904849, 0.533604905984294]);
    d = P([0.5490673726445683, 0.15991109487875477]);
    @test degree(gcd(a*d,b*d)) == 0
    @test degree(gcd(a*d, b*d, atol=sqrt(eps()))) > 0
    @test degree(gcd(a*d,b*d, method=:noda_sasaki)) == degree(d)
    @test_skip degree(gcd(a*d,b*d, method=:numerical)) == degree(d) # issues on some architectures (had test_skip)
    l,m,n = (4,4,4) #(5,5,5) # sensitive to choice of `rtol` in ngcd
    u,v,w = fromroots.(rand.((l,m,n)))
    @test degree(gcd(u*v, u*w, method=:numerical)) == degree(u)

    # Example of Zeng
    x = variable(P{Float64})
    p = (x+10)*(x^9 + x^8/3 + 1)
    q = (x+10)*(x^9 + x^8/7 - 6//7)

    @test degree(gcd(p,q)) == 0
    (@test degree(gcd(p,q, method=:noda_sasaki)) == 1)
    @test degree(gcd(p,q, method=:numerical)) == 1

    # more bits don't help Euclidean
    x = variable(P{BigFloat})
    p = (x+10)*(x^9 + x^8/3 + 1)
    q = (x+10)*(x^9 + x^8/7 - 6//7)
    @test degree(gcd(p,q)) == 0

    # Test 1 of Zeng
    x =  variable(P{Float64})
    alpha(j,n) = cos(j*pi/n)
    beta(j,n) = sin(j*pi/n)
    r1, r2 = 1/2, 3/2
    U(n) = prod( (x-r1*alpha(j,n))^2 + r1^2*beta(j,n)^2 for j in 1:n)
    V(n) = prod( (x-r2*alpha(j,n))^2 + r2^2*beta(j,n)^2 for j in 1:n)
    W(n) = prod( (x-r1*alpha(j,n))^2 + r1^2*beta(j,n)^2 for j in (n+1):2n)
    @testset for n in 2:2:20
        p = U(n) * V(n); q = U(n) * W(n)
        @test degree(gcd(p,q;  method=:numerical)) == degree(U(n))
    end

    # Test 5 of Zeng
    x =  variable(P{Float64})
    @testset for ms in ((2,1,1,0), (3,2,1,0), (4,3,2,1), (5,3,2,1), (9,6,4,2),
               (20, 14, 10, 5), (80,60,40,20), (100,60,40,20)
               )

        p = prod((x-i)^j for (i,j) in enumerate(ms))
        dp = derivative(p)
        @test degree(gcd(p,dp; method=:numerical)) == sum(max.(ms .- 1, 0))
    end

    # fussy pair
    x =  variable(P{Float64})
    @testset for n in (2,5,10,20,25,50, 100)
        p = (x-1)^n * (x-2)^n * (x-3)
        q = (x-1) * (x-2) * (x-4)
        a = Polynomials.ngcd(p, q)
        a.κ < 100 && @test degree(a.u) == 2
    end

    # check for fixed k
    p = fromroots(P, [2,3,4])
    q = fromroots(P, [3,4,5])
    out = Polynomials.ngcd(p,q)
    out1 = Polynomials.ngcd(p,q,1)
    out3 = Polynomials.ngcd(p,q,3)
    @test out.Θ <= out1.Θ
    @test out.Θ <= out3.Θ

    # check for correct output if degree p < degree q
    x = variable(P{Float64})
    p = -18.0 - 37.0*x - 54.0*x^2 - 36.0*x^3 - 16.0*x^4
    q = 2.0 + 5.0*x + 8.0*x^2 + 7.0*x^3 + 4.0*x^4 + 1.0*x^5
    out = Polynomials.ngcd(p,q)
    @test out.u * out.v ≈ p

    # check for canceling of x^k terms
    x = variable(P{Float64})
    p,q = x^2 + 1, x^2 - 1
    @testset for j ∈ 0:2
        @testset for k ∈ 0:j
            out = Polynomials.ngcd(x^j*p, x^k*q)
            @test out.u == x^k
        end
    end
end

@testset "Showing" begin

    # customized printing with printpoly
    function printpoly_to_string(args...; kwargs...)
        buf = IOBuffer()
        printpoly(buf, args...; kwargs...)
        return String(take!(buf))
    end

    p = Polynomial{Rational}([1, 4])
    @test sprint(show, p) == "Polynomial(1//1 + 4//1*x)"

    p = Polynomial{Rational{Int}}([1, 4])
    @test sprint(show, p) == "Polynomial(1//1 + 4//1*x)"

    @testset for P in (Polynomial, ImmutablePolynomial)
        p = P([1, 2, 3])
        @test sprint(show, p) == "$P(1 + 2*x + 3*x^2)"

        p = P([1.0, 2.0, 3.0])
        @test sprint(show, p) == "$P(1.0 + 2.0*x + 3.0*x^2)"

        p = P([1 + 1im, -2im])
        @test sprint(show, p) == "$P((1 + im) - 2im*x)"


        p = P([1,2,3,1])  # leading coefficient of 1
        @test repr(p) == "$P(1 + 2*x + 3*x^2 + x^3)"
        p = P([1.0, 2.0, 3.0, 1.0])
        @test repr(p) == "$P(1.0 + 2.0*x + 3.0*x^2 + 1.0*x^3)"
        p = P([1, im])
        @test repr(p) == "$P(1 + im*x)"
        p = P([1 + im, 1 - im, -1 + im, -1 - im])# minus signs
        @test repr(p) == "$P((1 + im) + (1 - im)x - (1 - im)x^2 - (1 + im)x^3)"
        p = P([1.0, 0 + NaN * im, NaN, Inf, 0 - Inf * im]) # handle NaN or Inf appropriately
        @test repr(p) == "$P(1.0 + NaN*im*x + NaN*x^2 + Inf*x^3 - Inf*im*x^4)"

        p = P([1,2,3])

        @test repr("text/latex", p) == "\$1 + 2\\cdot x + 3\\cdot x^{2}\$"
        p = P([1 // 2, 2 // 3, 1])
        @test repr("text/latex", p) == "\$\\frac{1}{2} + \\frac{2}{3}\\cdot x + x^{2}\$"
        p = P([complex(1,1),complex(0,1),complex(1,0),complex(1,1)])
        @test repr("text/latex", p) == "\$(1 + i) + i\\cdot x + x^{2} + (1 + i)x^{3}\$"

        @test printpoly_to_string(P([1,2,3], "y")) == "1 + 2*y + 3*y^2"
        @test printpoly_to_string(P([1,2,3], "y"), descending_powers = true) == "3*y^2 + 2*y + 1"
        @test printpoly_to_string(P([2, 3, 1], :z), descending_powers = true, offset = -2) == "1 + 3*z^-1 + 2*z^-2"
        @test printpoly_to_string(P([-1, 0, 1], :z), offset = -1, descending_powers = true) == "z - z^-1"
        @test printpoly_to_string(P([complex(1,1),complex(1,-1)]),MIME"text/latex"()) == "(1 + i) + (1 - i)x"
    end

    ## closed issues
    ## issue 275 with compact mult symbol
    p = Polynomial([1.234567890, 2.34567890])
    io=IOBuffer(); printpoly(io, p, compact=true); @test String(take!(io)) == "1.23457 + 2.34568*x"
    io=IOBuffer(); printpoly(io, p, compact=true, mulsymbol=""); @test String(take!(io)) == "1.23457 + 2.34568x"

    ## issue 278 with complex
    @test printpoly_to_string(Polynomial([1 + im, 1, 2, im, 2im, 1+im, 1-im])) == "(1 + im) + x + 2*x^2 + im*x^3 + 2im*x^4 + (1 + im)x^5 + (1 - im)x^6"

    ## issue #320 (fix was broken)
    @test printpoly_to_string(Polynomial(BigInt[1,0,1], :y)) == "1 + y^2"

    # negative indices
    @test printpoly_to_string(LaurentPolynomial([-1:3;], -2)) == "-x⁻² + 1 + 2*x + 3*x²"
    @test printpoly_to_string(SparsePolynomial(Dict(.=>(-2:2, -1:3)))) == "-x^-2 + 1 + 2*x + 3*x^2"
end

@testset "Plotting" begin
    p = fromroots([-1, 1]) # x^2 - 1
    rec = apply_recipe(Dict{Symbol,Any}(), p)
    @test rec[1].plotattributes[:label] == "-1 + x^2"
    @test rec[1].plotattributes[:xlims] == (-1.4, 1.4)


    rec = apply_recipe(Dict{Symbol,Any}(), p, -1, 1)
    @test rec[1].plotattributes[:label] == "-1 + x^2"

    p = ChebyshevT([1,1,1])
    rec = apply_recipe(Dict{Symbol,Any}(), p)
    @test !isnothing(match(r"T_0", rec[1].plotattributes[:label]))
    @test rec[1].plotattributes[:xlims] == (-1.0, 1.0) # uses domain(p)
end

@testset "Promotion"  begin

    # Test different types work together
    @testset for P₁ in Ps
        @testset for   P₂ in Ps
            p₁, p₂ = P₁(rand(1:5, 4)), P₂(rand(1:5, 5))
            p₁ + p₂
            p₁ * p₂

            p₁, p₂ = P₁(rand(1:5, 4)), P₂(5) # constant
            p₁ + p₂
            p₁ * p₂

            p₁, p₂ = P₁(rand(1:5, 4)), P₂(5, :y) # constant, but wrong variable
            if !(promote_type(P₁, P₂) <: Polynomial || promote_type(P₁, P₂) <: Polynomials.StandardBasisPolynomial)
                p₁ + p₂
                p₁ * p₂
            end
        end
    end

    # P{T}(vector{S}) will promote to P{T}
    @testset for Ts in ((Int32, Int,  BigInt),
               (Int,  Rational{Int}, Float64),
               (Float32, Float64, BigFloat)
              )

        n = length(Ts)
        @testset for i in 1:n-1
            T1,T2 = Ts[i],Ts[i+1]
            @testset for P in Ps
                P <: FactoredPolynomial && continue
                if !(P ∈ (ImmutablePolynomial, ImmutableDensePolynomial{StandardBasis}))
                    p = P{T2}(T1.(rand(1:3,3)))
                    @test typeof(p) == P{T2, :x}
                else
                    N = 3
                    p = P{T2}(T1.(rand(1:3,N)))
                    @test typeof(p) == P{T2,:x, N}
                end
            end

        end
    end

    # test P{T}(...) is P{T} (not always the case for FactoredPolynomial)
    @testset for P in Ps
        P <: FactoredPolynomial && continue
        if !(P ∈ (ImmutablePolynomial, ImmutableDensePolynomial{StandardBasis}))
            @testset for T in (Int32, Int64, BigInt)
                p₁ =  P{T}(Float64.(rand(1:3,5)))
                @test typeof(p₁) == P{T,:x} # conversion works
                @test_throws InexactError  P{T}(rand(5))
            end
        else
            @testset for T in (Int32, Int64, BigInt)
                N = 5
                p₁ =  P{T}(Float64.(rand(1:3,5)))
                @test typeof(p₁) == P{T,:x,5} # conversion works
                @test_throws InexactError  P{T}(rand(5))
            end
        end
    end
end


@testset "SparsePolynomial" begin
    @test Polynomials.minimumexponent(SparsePolynomial) == typemin(Int)
    @test Polynomials.minimumexponent(SparsePolynomial{Float64}) == typemin(Int)
    @test Polynomials.minimumexponent(SparsePolynomial{Float64, :y}) == typemin(Int)
    @testset "empty" begin
        p = SparsePolynomial(Float64[0])
        @test eltype(p) == Float64
        @test eltype(keys(p)) == Int
        @test eltype(values(p)) == Float64
        @test collect(p) == Float64[]
        @test collect(keys(p)) == Int[]
        @test collect(values(p)) == Float64[]
        @test p == Polynomial(0)
    end
    @testset "negative indices" begin
        d = Dict(-2=>4, 5=>10)
        p = SparsePolynomial(d)
        @test length(p) == 8
        @test firstindex(p) == -2
        @test lastindex(p) == 5
        @test eachindex(p) == -2:5
        q = LaurentPolynomial(p)
        @test p == q
        @test SparsePolynomial(q) == p
        @test_throws ArgumentError Polynomial(p)
    end
end

@testset "LaurentPolynomial" begin
    @test Polynomials.minimumexponent(LaurentPolynomial) == typemin(Int)
    @test Polynomials.minimumexponent(LaurentPolynomial{Float64}) == typemin(Int)
    @test Polynomials.minimumexponent(LaurentPolynomial{Float64, :y}) == typemin(Int)
end


# Chain rules
using ChainRulesTestUtils

@testset "Test frule and rrule" begin
    p = Polynomial([1,2,3,4])
    dp = derivative(p)

    test_scalar(p, 1.0; check_inferred=true)
end
