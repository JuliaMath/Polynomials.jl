# assert file to test polynomial implementation
using Test
using LinearAlgebra
using Polynomials
using SpecialFunctions

import SparseArrays: sparse, nnz

@testset "Polynomial" begin include("Polynomial.jl") end
@testset "Old" begin include("old.jl") end
