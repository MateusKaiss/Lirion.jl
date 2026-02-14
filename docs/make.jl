import Documenter
using Lirion

Documenter.makedocs(
    sitename = "Lirion.jl",
    authors = "Mateus Kaiss and contributors.",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        prettyurls = false,
        sidebar_sitename = false,
        assets = ["assets/custom.css"], # se você criou o custom.css
    ),
    modules = [Lirion],
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
    ],
)

Documenter.deploydocs(
    repo = "github.com/MateusKaiss/Lirion.jl.git",
)
