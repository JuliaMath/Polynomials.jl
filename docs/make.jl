using Polynomials, Documenter

makedocs(modules = [Polynomials],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Polynomials.jl",
    authors = "Jameson Nash, Keno Fischer, and other contributors",
    pages = [
        "Manual" => "index.md",
    ],
)

deploydocs(repo = "github.com/JuliaMath/Polynomials.jl.git")
