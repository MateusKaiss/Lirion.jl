import Documenter
using Lirion

Documenter.makedocs(
    sitename = "Lirion.jl",
    authors = "Mateus Kaiss and contributors.",
    format = Documenter.HTML(;
        mathengine = Documenter.MathJax(),
        prettyurls = false,
        sidebar_sitename = false,
        collapselevel = 1,
    ),
    modules = [Lirion],
    checkdocs = :none,
    pages = [
        "Introduction" => [
            "Introduction" => "index.md","installation.md","quick_start.md"
        ],
        "Performance" => "performance.md",
    ],
)

Documenter.deploydocs(
    repo = "github.com/MateusKaiss/Lirion.jl.git",
)
