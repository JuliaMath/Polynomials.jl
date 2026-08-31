using Documenter
using DocumenterCodeBlocks, DocumenterCitations, DocumenterLandingPage
using Polynomials


ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"


DocMeta.setdocmeta!(Polynomials, :DocTestSetup, :(using Polynomials); recursive = true)

makedocs(
    modules = [Polynomials],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Polynomials.jl",
    authors = "Jameson Nash, Keno Fischer, Miles Lucas, John Verzani, and other contributors",
    pages = [
        "Home" => "index.md",
        "Reference/API" => "reference.md",
        "Polynomial Types" => [
            "Polynomial" => "polynomials/polynomial.md",
            "Chebyshev" => "polynomials/chebyshev.md",
        ],
        "Extending" => "extending.md",
    ],
    plugins = [LandingPage(),
               CodeBlocks(),
               CitationBibliography("src/refs.bib")],
    doctestfilters = [
        r"(?<=\d\.\d{12})\d+", # Ignore any digit after the 12th decimal place
    ],
    warnonly = [:cross_references, :missing_docs],
    checkdocs=:exports,
)

deploydocs(repo = "github.com/JuliaMath/Polynomials.jl.git")
