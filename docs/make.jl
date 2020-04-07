using Documenter
using Polynomials

DocMeta.setdocmeta!(Polynomials, :DocTestSetup, :(using Polynomials); recursive = true)

makedocs(
    modules = [Polynomials],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Polynomials.jl",
    authors = "Jameson Nash, Keno Fischer, and other contributors",
    pages = [
        "Home" => "index.md",
        "Reference/API" => "reference.md",
        "Polynomial Types" => [
            "Polynomial" => "polynomials/polynomial.md",
            "Chebyshev" => "polynomials/chebyshev.md",
        ],
        "Extending" => "extending.md",
    ],
)

deploydocs(repo = "github.com/JuliaMath/Polynomials.jl.git")
