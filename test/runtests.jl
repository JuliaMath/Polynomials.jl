# assert file to test polynomial implementation
using Test
using LinearAlgebra
using Polynomials
using SpecialFunctions
using RecipesBase: apply_recipe

import SparseArrays: sparse, nnz
using OffsetArrays

@testset "Standard basis" begin include("StandardBasis.jl") end
@testset "ChebyshevT" begin include("ChebyshevT.jl") end
@testset "Rational functions" begin include("rational-functions.jl") end
@testset "Poly, Pade (compatability)" begin include("Poly.jl") end
