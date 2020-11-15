using Documenter, DiffEqBayes

makedocs(
    sitename="DiffEqBayes.jl",
    authors="Chris Rackauckas et al.",
    clean=true,
    doctest=false,
    modules=[DiffEqBayes],

    format=Documenter.HTML(assets=["assets/favicon.ico"],
                           canonical="https://diffeqbayes.sciml.ai/stable/"),

    pages=[
        "DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations" => "index.md",
        "Methods" => "methods.md",
        "Examples" => "examples.md"
    ]
)

deploydocs(
    repo="github.com/SciML/DiffEqBayes.jl.git";
    push_preview=true
)
