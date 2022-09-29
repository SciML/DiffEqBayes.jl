"""
$(DocStringExtensions.README)
"""
module DiffEqBayes

# Import dependencies alphabetically
using DiffEqBase
using Distances
using Distributions
using Distributions
using DocStringExtensions
using DocStringExtensions
using LinearAlgebra
using MacroTools
using ModelingToolkit
using Optim
using Parameters
using Random
using RecursiveArrayTools
using Requires
using StanSample
using Turing

STANDARD_PROB_GENERATOR(prob, p) = remake(prob; u0 = eltype(p).(prob.u0), p = p)
function STANDARD_PROB_GENERATOR(prob::EnsembleProblem, p)
    EnsembleProblem(remake(prob.prob; u0 = eltype(p).(prob.prob.u0), p = p))
end

include("turing_inference.jl")
# include("abc_inference.jl")
include("stan_string.jl")
include("stan_inference.jl")

function __init__()
    @require DynamicHMC="bbc10e6e-7c05-544b-b16e-64fede858acb" begin
        using .DynamicHMC, TransformVariables, LogDensityProblems, TransformedLogDensities
        include("dynamichmc_inference.jl")
        export dynamichmc_inference
    end
end

export turing_inference, stan_inference, abc_inference
end # module
