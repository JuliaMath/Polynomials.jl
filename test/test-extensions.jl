using FFTW
using Makie
using ChainRulesCore
using RecipesBase: apply_recipe

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
