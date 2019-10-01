using Polynomials, Documenter

makedocs(modules = [Polynomials],
    clean = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Polynomials.jl",
    authors = "Jameson Nash, Keno Fischer, and other contributors",
    pages = [
        "Manual" => "index.md",
    ],
)

deploydocs("github.com/JuliaMath/Polynomials.jl.git")
