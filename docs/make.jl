using Documenter, DiffEqBayes

makedocs(sitename="DiffEqBayes.jl",
    modules=[DiffEqBayes],
    clean=true,doctest=false,
    format = Documenter.HTML(#analytics = "UA-90474609-3",
    assets = ["assets/favicon.ico"])
)

deploydocs(
   repo = "github.com/SciML/DiffEqBayes.jl.git";
   push_preview = true
)