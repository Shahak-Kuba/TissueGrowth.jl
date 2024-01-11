using Documenter
using TissueGrowth

push!(LOAD_PATH,"../src/")

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
