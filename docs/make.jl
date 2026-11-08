# =============================================================================
# Lirion.jl — Documentation 
# =============================================================================

import Documenter
using Documenter
using Lirion

bib = CitationBibliography(
    joinpath(@__DIR__, "refs.bib");
    style = :numeric,          # :numeric | :authoryear | :alpha
)

Documenter.makedocs(
    sitename = "Lirion.jl",
    authors  = "Mateus Kaiss and contributors.",
    format = Documenter.HTML(;
        mathengine = Documenter.MathJax3(),
        prettyurls = false,
        sidebar_sitename = false,
        collapselevel = 1,
    ),
    modules   = [Lirion],
    checkdocs = :none,   
    plugins   = [bib],       
    pages = [
        "Introduction" => "index.md",
        "Getting Started" => [
            "installation.md",
            "Exact formulations" => [
                "AC Polar" => "usage/exact/ac_polar.md",
                "ACR"      => "usage/exact/acr.md",
                "IVR"      => "usage/exact/ivr.md",
            ],
            "Relaxations" => [
                "SOCWR"  => "usage/relaxations/socwr.md",
                "QC"     => "usage/relaxations/qc.md",
                "SDPWRM" => "usage/relaxations/sdpwrm.md",
            ],
            "Approximations" => [
                "DC"   => "usage/approximations/dc.md",
                "LPAC" => "usage/approximations/lpac.md",
                "NFA"  => "usage/approximations/nfa.md",
                "BFA"  => "usage/approximations/bfa.md",
            ],
            "Tunable settings (centralized)" => "usage/tunable_settings.md",
            "Distributed" => [
                "Running a distributed solve" => "usage/distributed/overview.md",
                "Coupling and sectioning"     => "usage/distributed/coupling_sectioning.md",
                "ADMM settings"               => "usage/distributed/admm.md",
                "ALADIN settings"             => "usage/distributed/aladin.md",
            ],
            "Results" => "usage/results.md",
        ],
        "Theory" => [
            "Exact formulations" => [
                "AC Polar"              => "teoria/exatas/ac_polar.md",
                "AC Rectangular (ACR)"  => "teoria/exatas/acr.md",
                "IVR"                   => "teoria/exatas/ivr.md",
            ],
            "Relaxations" => [
                "SOCWR (Jabr)" => "teoria/relaxacoes/socwr.md",
                "QC"           => "teoria/relaxacoes/qc.md",
                "SDPWRM"       => "teoria/relaxacoes/sdpwrm.md",
            ],
            "Approximations" => [
                "DC"                => "teoria/aproximacoes/dc.md",
                "LPAC"              => "teoria/aproximacoes/lpac.md",
                "NFA"               => "teoria/aproximacoes/nfa.md",
                "BFA (LinDistFlow)" => "teoria/aproximacoes/bfa.md",
            ],
            "Distributed" => [
                "Sectioning and coupling" => "teoria/distribuido/seccionamento_acoplamento.md",
                "ADMM"   => "teoria/distribuido/admm.md",
                "ALADIN" => "teoria/distribuido/aladin.md",
            ],
        ],
        "Performance" => "performance.md",
        "References" => "references.md",
    ],
)

# --- Deployment: INTENTIONALLY DISABLED in the scaffold ----------------------
# Do NOT uncomment here. Running this would push the built site to the
# gh-pages branch of origin. Re-enable only in the release, after confirming
# the repo slug and that a deploy key / GITHUB_TOKEN is configured.
#
# Documenter.deploydocs(
#     repo         = "github.com/MateusKaiss/Lirion.jl.git",
#     push_preview = false,
# )
