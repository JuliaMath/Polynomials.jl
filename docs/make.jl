using Polynomials, Documenter

makedocs(
    modules = [Polynomials],
    clean = false,
    format = :html,
    sitename = "Polynomials.jl",
    authors = "Jameson Nash, Keno Fischer, and other contributors",
    pages = [
        "Manual" => "index.md",
    ],
)

deploydocs(
    julia = "nightly",
    repo = "github.com/JuliaMath/Polynomials.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
