using Documenter
using TissueGrowth

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "TissueGrowth.jl",
    modules = [TissueGrowth],
    pages = [
        "Home" => "index.md",
        "ForUser.md",
        "Equations.md",
        "Example.md",
        "PDEsolver.md"
    ],
warnonly = Documenter.except(),
format = Documenter.HTML(prettyurls = false)
)

deploydocs(;
    repo = "github.com/Shahak-Kuba/TissueGrowth.jl",
    devbranch="main",
)
