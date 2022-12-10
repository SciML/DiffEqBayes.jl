"""
$(DocStringExtensions.README)
"""
module DiffEqBayes

using DiffEqBase, Distributions, Turing, MacroTools
using RecursiveArrayTools, ModelingToolkit, LinearAlgebra
using Parameters, Distributions, Optim, Requires
using Distances, DocStringExtensions, Random, StanSample
using DynamicHMC, TransformVariables, LogDensityProblems, TransformedLogDensities

STANDARD_PROB_GENERATOR(prob, p) = remake(prob; u0 = eltype(p).(prob.u0), p = p)
function STANDARD_PROB_GENERATOR(prob::EnsembleProblem, p)
    EnsembleProblem(remake(prob.prob; u0 = eltype(p).(prob.prob.u0), p = p))
end

include("turing_inference.jl")
# include("abc_inference.jl")
include("stan_string.jl")
include("stan_inference.jl")
include("dynamichmc_inference.jl")

export turing_inference, stan_inference, dynamichmc_inference
end # module
