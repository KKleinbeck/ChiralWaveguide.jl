using Documenter, ChiralWaveguide

makedocs(
  sitename = "ChiralWaveguide",
  modules = [ChiralWaveguide],
  format = Documenter.HTML(assets = ["assets/favicon.ico"]),
  pages = [
    "index.md",
    # "theory.md",
    "Examples" => [
      "examples/SinglePhotonScattering.md",
      "examples/CustomSystems.md"
    ],
    "api.md"
  ]
)

deploydocs(
  repo = "github.com/KKleinbeck/ChiralWaveguide.jl.git",
)
