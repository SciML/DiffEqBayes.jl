using Documenter, DiffEqBayes

makedocs(
    sitename="DiffEqBayes.jl",
    authors="#",
    clean=true,
    doctest=false,
    modules=[DiffEqBayes],

    format=Documenter.HTML(assets=["assets/favicon.ico"],
                           canonical="#"),

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
