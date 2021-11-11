using Documenter
using EvolutionaryModelingTools

makedocs(;
    sitename="EvolutionaryModelingTools.jl",
    pages=["index.md", "example.md", "reference.md"],
)

deploydocs(;
    repo="github.com/wangl-cc/EvolutionaryModelingTools.jl.git",
    push_preview=true,
)
