var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Polynomials.jl-1",
    "page": "Manual",
    "title": "Polynomials.jl",
    "category": "section",
    "text": "Polynomials.jl is a Julia package that provides basic arithmetic, integration, differentiation, evaluation, and root finding over dense univariate polynomials.To install the package, runPkg.add(\"Polynomials\")The package can then be loaded into the current session usingusing PolynomialsDocTestSetup = :( using Polynomials )"
},

{
    "location": "index.html#Polynomials.Poly",
    "page": "Manual",
    "title": "Polynomials.Poly",
    "category": "Type",
    "text": "Poly{T<:Number}(a::AbstractVector{T}, [x])\n\nConstruct a polynomial from its coefficients a, lowest order first, optionally in terms of the given variable x. x can be a character, symbol, or string.\n\nIf p = a_n x^n + ldots + a_2 x^2 + a_1 x + a_0, we construct this through Poly([a_0, a_1, ..., a_n]).\n\nThe usual arithmetic operators are overloaded to work with polynomials as well as with combinations of polynomials and scalars. However, operations involving two polynomials of different variables causes an error.\n\nExamples\n\njulia> Poly([1, 0, 3, 4])\nPoly(1 + 3⋅x^2 + 4⋅x^3)\n\njulia> Poly([1, 2, 3], :s)\nPoly(1 + 2⋅s + 3⋅s^2)\n\njulia> a = Poly([1, 2, 3], :x); b = Poly([1, 2, 3], :s);\n\njulia> a + b\nERROR: Polynomials must have same variable\n...\n\njulia> p = Poly([1, 2])\nPoly(1 + 2⋅x)\n\njulia> q = Poly([1, 0, -1])\nPoly(1 - x^2)\n\njulia> 2p\nPoly(2 + 4⋅x)\n\njulia> 2 + p\nPoly(3 + 2⋅x)\n\njulia> p - q\nPoly(2⋅x + x^2)\n\njulia> p * q\nPoly(1 + 2⋅x - x^2 - 2⋅x^3)\n\njulia> q / 2\nPoly(0.5 - 0.5⋅x^2)\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.poly",
    "page": "Manual",
    "title": "Polynomials.poly",
    "category": "Function",
    "text": "poly(r)\n\nConstruct a polynomial from its roots. Compare this to the Poly type constructor, which constructs a polynomial from its coefficients.\n\nIf r is a vector, the constructed polynomial is (x - r_1) (x - r_2) cdots (x - r_n). If r is a matrix, the constructed polynomial is (x - e_1) cdots (x - e_n), where e_i is the ith eigenvalue of r.\n\nExamples\n\njulia> poly([1, 2, 3])   # The polynomial (x - 1)(x - 2)(x - 3)\nPoly(-6 + 11⋅x - 6⋅x^2 + x^3)\n\njulia> poly([1 2; 3 4])  # The polynomial (x - 5.37228)(x + 0.37228)\nPoly(-1.9999999999999998 - 5.0⋅x + 1.0⋅x^2)\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.degree",
    "page": "Manual",
    "title": "Polynomials.degree",
    "category": "Function",
    "text": "degree(p::Poly)\n\nReturn the degree of the polynomial p, i.e. the highest exponent in the polynomial that has a nonzero coefficient.\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.coeffs",
    "page": "Manual",
    "title": "Polynomials.coeffs",
    "category": "Function",
    "text": "coeffs(p::Poly)\n\nReturn the coefficient vector [a_0, a_1, ..., a_n] of a polynomial p.\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.variable",
    "page": "Manual",
    "title": "Polynomials.variable",
    "category": "Function",
    "text": "variable(p::Poly)\nvariable([T::Type,] var)\nvariable()\n\nReturn the indeterminate of a polynomial, i.e. its variable, as a Poly object. When passed no arguments, this is equivalent to variable(Float64, :x).\n\nExamples\n\njulia> variable(Poly([1, 2], :x))\nPoly(x)\n\njulia> variable(:y)\nPoly(1.0⋅y)\n\njulia> variable()\nPoly(1.0⋅x)\n\njulia> variable(Float32, :x)\nPoly(1.0f0⋅x)\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.printpoly",
    "page": "Manual",
    "title": "Polynomials.printpoly",
    "category": "Function",
    "text": "printpoly(io::IO, p::Poly, mimetype = MIME\"text/plain\"(); descending_powers=false)\n\nPrint a human-readable representation of the polynomial p to io. The MIME types \"text/plain\" (default), \"text/latex\", and \"text/html\" are supported. By default, the terms are in order of ascending powers, matching the order in coeffs(p); specifying descending_powers=true reverses the order.\n\nExamples\n\njulia> printpoly(STDOUT, Poly([1,2,3], :y))\n1 + 2*y + 3*y^2\njulia> printpoly(STDOUT, Poly([1,2,3], :y), descending_powers=true)\n3*y^2 + 2*y + 1\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.polyval",
    "page": "Manual",
    "title": "Polynomials.polyval",
    "category": "Function",
    "text": "polyval(p::Poly, x::Number)\n\nEvaluate the polynomial p at x using Horner's method. Poly objects are callable, using this function.\n\nExamples\n\njulia> p = Poly([1, 0, -1])\nPoly(1 - x^2)\n\njulia> polyval(p, 1)\n0\n\njulia> p(1)\n0\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.polyint",
    "page": "Manual",
    "title": "Polynomials.polyint",
    "category": "Function",
    "text": "polyint(p::Poly, k::Number=0)\n\nIntegrate the polynomial p term by term, optionally adding a constant term k. The order of the resulting polynomial is one higher than the order of p.\n\nExamples\n\njulia> polyint(Poly([1, 0, -1]))\nPoly(1.0⋅x - 0.3333333333333333⋅x^3)\n\njulia> polyint(Poly([1, 0, -1]), 2)\nPoly(2.0 + 1.0⋅x - 0.3333333333333333⋅x^3)\n\n\n\n\n\npolyint(p::Poly, a::Number, b::Number)\n\nCompute the definite integral of the polynomial p over the interval [a,b].\n\nExamples\n\njulia> polyint(Poly([1, 0, -1]), 0, 1)\n0.6666666666666667\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.polyder",
    "page": "Manual",
    "title": "Polynomials.polyder",
    "category": "Function",
    "text": "polyder(p::Poly, k=1)\n\nCompute the kth derivative of the polynomial p.\n\nExamples\n\njulia> polyder(Poly([1, 3, -1]))\nPoly(3 - 2⋅x)\n\njulia> polyder(Poly([1, 3, -1]), 2)\nPoly(-2)\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.polyfit",
    "page": "Manual",
    "title": "Polynomials.polyfit",
    "category": "Function",
    "text": "polyfit(x, y, n=length(x)-1, sym=:x)\n\nFit a polynomial of degree n through the points specified by x and y, where n <= length(x) - 1, using least squares fit. When n=length(x)-1 (the default), the interpolating polynomial is returned. The optional fourth argument can be used to specify the symbol for the returned polynomial.\n\nExamples\n\njulia> xs = linspace(0, pi, 5);\n\njulia> ys = map(sin, xs);\n\njulia> polyfit(xs, ys, 2)\nPoly(-0.004902082150108854 + 1.242031920509868⋅x - 0.39535103925413095⋅x^2)\n\n\n\n\n\n"
},

{
    "location": "index.html#Polynomials.roots",
    "page": "Manual",
    "title": "Polynomials.roots",
    "category": "Function",
    "text": "roots(p::Poly)\n\nReturn the roots (zeros) of p, with multiplicity. The number of roots returned is equal to the order of p. The returned roots may be real or complex.\n\nExamples\n\njulia> roots(Poly([1, 0, -1]))\n2-element Array{Float64,1}:\n -1.0\n  1.0\n\njulia> roots(Poly([1, 0, 1]))\n2-element Array{Complex{Float64},1}:\n 0.0+1.0im\n 0.0-1.0im\n\njulia> roots(Poly([0, 0, 1]))\n2-element Array{Float64,1}:\n 0.0\n 0.0\n\njulia> roots(poly([1,2,3,4]))\n4-element Array{Float64,1}:\n 4.0\n 3.0\n 2.0\n 1.0\n\n\n\n\n\n"
},

{
    "location": "index.html#Functions-1",
    "page": "Manual",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = PolynomialsPoly\npoly\ndegree\ncoeffs\nvariable\nprintpoly\npolyval\npolyint\npolyder\npolyfit\nroots"
},

]}
