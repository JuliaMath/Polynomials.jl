# assert file to test polynomial implementation
using Test
using LinearAlgebra
using Polynomials
using SpecialFunctions

import SparseArrays: sparse, nnz

@testset "Polynomial" begin include("Polynomial.jl") end
@testset "Deprecations" begin include("deprecated.jl") end
@testset "Poly (deprecaterd)" begin include("Poly.jl") end
