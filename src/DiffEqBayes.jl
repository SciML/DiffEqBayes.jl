"""
$(DocStringExtensions.README)
"""
module DiffEqBayes

using DocStringExtensions
using DiffEqBase, Distributions, Turing, MacroTools
using RecursiveArrayTools, ModelingToolkit, LinearAlgebra
using Parameters, Distributions, Optim, Requires
using Distances, DocStringExtensions, Random, StanSample

STANDARD_PROB_GENERATOR(prob,p) = remake(prob;u0=eltype(p).(prob.u0),p=p)
STANDARD_PROB_GENERATOR(prob::EnsembleProblem,p) = EnsembleProblem(remake(prob.prob;u0=eltype(p).(prob.prob.u0),p=p))

include("turing_inference.jl")
# include("abc_inference.jl")
include("stan_string.jl")
include("stan_inference.jl")

function __init__()
    @require DynamicHMC="bbc10e6e-7c05-544b-b16e-64fede858acb" begin
        using .DynamicHMC, TransformVariables, LogDensityProblems
        include("dynamichmc_inference.jl")
        export dynamichmc_inference
    end
end

export turing_inference, stan_inference ,abc_inference
end # module
