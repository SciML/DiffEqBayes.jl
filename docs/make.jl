using Documenter, DiffEqBayes

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"
include("pages.jl")

makedocs(sitename = "DiffEqBayes.jl",
         authors = "Chris Rackauckas, Vaibhav Kumar Dixit et al.",
         clean = true,
         doctest = false,
         modules = [DiffEqBayes],
         strict = [
             :doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # Other available options are
             # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
         ],
         format = Documenter.HTML(assets = ["assets/favicon.ico"],
                                  canonical = "https://docs.sciml.ai/DiffEqBayes/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/DiffEqBayes.jl.git";
           push_preview = true)
