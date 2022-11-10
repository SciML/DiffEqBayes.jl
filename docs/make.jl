using Documenter, DiffEqBayes

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "DiffEqBayes.jl",
         authors = "Chris Rackauckas, Vaibhav Kumar Dixit et al.",
         clean = true,
         doctest = false,
         modules = [DiffEqBayes],
         format = Documenter.HTML(assets = ["assets/favicon.ico"],
                                  canonical = "https://docs.sciml.ai/DiffEqBayes/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/DiffEqBayes.jl.git";
           push_preview = true)
