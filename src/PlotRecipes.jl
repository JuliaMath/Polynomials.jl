using RecipesBase

@recipe function poly_recipe(p::Poly, xs = polyinterval(p))
    buff = IOBuffer()
    printpoly(buff, p)
    label --> String(take!(buff))

    ys = map(x -> polyval(p, x), xs)
    xs, ys
end
