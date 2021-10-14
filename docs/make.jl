using Documenter, ChiralWaveguide

makedocs(
  sitename = "ChiralWaveguide",
  modules = [ChiralWaveguide],
  pages = [
    "index.md",
    # "theory.md",
    "Examples" => [
      "examples/SinglePhotonScattering.md"
    ],
    "api.md"
  ]
)

deploydocs(
  repo = "github.com/KKleinbeck/ChiralWaveguide.jl.git",
)
