# assert file to test polynomial implementation
using Test
using LinearAlgebra
using Polynomials
using SpecialFunctions
using RecipesBase: apply_recipe

import SparseArrays: sparse, nnz

@testset "Polynomial" begin include("Polynomial.jl") end
#@testset "Immutable" begin include("Immutable.jl") end
@testset "ChebyshevT" begin include("ChebyshevT.jl") end
@testset "Poly, Pade (compatability)" begin include("Poly.jl") end
