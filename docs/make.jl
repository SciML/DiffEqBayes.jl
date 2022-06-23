using Documenter, DiffEqBayes

include("pages.jl")

makedocs(sitename = "DiffEqBayes.jl",
         authors = "Chris Rackauckas, Vaibhav Kumar Dixit et al.",
         clean = true,
         doctest = false,
         modules = [DiffEqBayes],
         format = Documenter.HTML(assets = ["assets/favicon.ico"],
                                  canonical = "https://diffeqbayes.sciml.ai/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/DiffEqBayes.jl.git";
           push_preview = true)
