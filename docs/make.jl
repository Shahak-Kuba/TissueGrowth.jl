using Documenter
using TissueGrowth

makedocs(
    sitename = "TissueGrowth.jl",
    modules = [TissueGrowth],
    pages = [
        "Home" => "index.md"
    ],
warnonly = Documenter.except()
)

deploydocs(;
    repo = "github.com/Shahak-Kuba/TissueGrowth.jl",
    devbranch="main",
)