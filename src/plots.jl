using RecipesBase

function poly_interval(p::AbstractPolynomial)

    # use  restricted domain, if finite
    A,B =  first(domain(p)), last(domain(p))
    if !isinf(A) && !isinf(B)
        # if isopen(domain(p))
        #     Delta = (B-A)/100
        #     A += Delta
        #     B -= Delta
        # end
        return A:(B-A)/100:B
    end


    # Find points of interest
    zero_pts = roots(p)
    crit_pts = roots(derivative(p, 1))
    infl_pts = roots(derivative(p, 2))
    pts = sort([ real(pt) for pt in [zero_pts; crit_pts; infl_pts] if isreal(pt) ])
    # Choose a range that shows all interesting points with some margin
    min_x, max_x = length(pts) > 0 ? (pts[1], pts[end]) : (-1, 1)
    d = max(max_x - min_x, 1)
    a = min_x - d / 5
    b = max_x + d / 5

    Delta  = b -  a

    return a:Delta/50:b
end

poly_label(p::AbstractPolynomial) = sprint(printpoly, p)

@recipe function f(p::AbstractPolynomial, x = poly_interval(p))
    label --> poly_label(p)
    x, p.(x)
end

@recipe function f(p::AbstractPolynomial, a, b)
    label --> poly_label(p)
    step = (b - a) / 100
    xs = a:step:b
    ys = p.(xs)
    xs, ys
end
