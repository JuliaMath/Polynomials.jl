# assert file to test polynomial implementation
using Test
using LinearAlgebra
using Polynomials
using SpecialFunctions
using RecipesBase: apply_recipe
using OffsetArrays

import SparseArrays: sparse, nnz

@testset "Standard basis" begin include("StandardBasis.jl") end
@testset "ChebyshevT" begin include("ChebyshevT.jl") end
@testset "Poly, Pade (compatability)" begin include("Poly.jl") end
