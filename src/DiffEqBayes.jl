module DiffEqBayes
using DiffEqBase, Mamba, Stan, Distributions, Turing, MacroTools
using Compat

include("stan_inference.jl")
include("turing_inference.jl")

export StanModel, stan_inference, turing_inference

end # module
