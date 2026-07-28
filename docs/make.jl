using Documenter
using WeberElectrodynamics

DocMeta.setdocmeta!(
    WeberElectrodynamics,
    :DocTestSetup,
    :(using WeberElectrodynamics);
    recursive = true,
)

makedocs(
    sitename = "WeberElectrodynamics.jl",
    authors = "Samer",
    modules = [WeberElectrodynamics],
    repo = "github.com/WeberElectrodynamics/WeberElectrodynamics.jl.git",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://WeberElectrodynamics.github.io/WeberElectrodynamics.jl",
        assets = String[],
        repolink = "https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quickstart.md",
        "The Weber Hamiltonian" => "hamiltonian.md",
        "API Reference" => [
            "System" => "api/system.md",
            "Problem" => "api/problem.md",
            "Solver" => "api/solver.md",
            "Callbacks" => "api/callbacks.md",
            "Statistics" => "api/statistics.md",
            "Visualization" => "api/visualization.md",
        ],
        "Internals" => "internals.md",
        "Advanced" => [
            "Regularization" => "regularization.md",
            "Custom Hamiltonians" => "custom_hamiltonians.md",
        ],
        "Theory" => "theory.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
)

if get(ENV, "GITHUB_EVENT_NAME", "") != "pull_request"
    deploydocs(
        repo = "github.com/WeberElectrodynamics/WeberElectrodynamics.jl.git",
        devbranch = "main",
        push_preview = false,
    )
end
