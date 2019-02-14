using RecipesBase


function poly_interval(p::Poly)
    # Find points of interest
    zero_pts = roots(p)
    crit_pts = roots(polyder(p, 1))
    infl_pts = roots(polyder(p, 2))
    pts = sort([ real(pt) for pt in [zero_pts; crit_pts; infl_pts] if isreal(pt) ])

    # Choose a range that shows all interesting points with some margin
    min_x, max_x = length(pts) > 0 ? (pts[1], pts[end]) : (-1, 1)
    diff = max(max_x - min_x, 1)
    a = min_x - diff/2
    b = max_x + diff/2

    return a : diff/50 : b
end


function poly_label(p::Poly)
    buff = IOBuffer()
    printpoly(buff, p)
    String(take!(buff))
end


@recipe function poly_recipe(p::Poly, range = poly_interval(p))
    label --> poly_label(p)
    range, map(x -> polyval(p, x), range)
end


@recipe function poly_recipe(p::Poly, a, b)
    label --> poly_label(p)
    step = (b-a)/100
    xs = a:step:b
    ys = map(x -> polyval(p, x), xs)
    xs, ys
end
