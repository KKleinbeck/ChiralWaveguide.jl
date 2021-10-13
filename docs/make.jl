using Documenter, ChiralWaveguide

makedocs(
    sitename = "ChiralWaveguide",
    pages = ["index.md", "api.md"]
)

deploydocs(
    repo = "github.com/KKleinbeck/ChiralWaveguide.jl.git",
)
