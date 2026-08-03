"""
$(DocStringExtensions.README)
"""
module DiffEqBayes

using DiffEqBase: DiffEqBase, EnsembleProblem, SciMLBase, remake
using Distances: Distances
using DocStringExtensions: DocStringExtensions, FIELDS, SIGNATURES, TYPEDEF
using DynamicHMC: DynamicHMC, mcmc_with_warmup
using LinearAlgebra: LinearAlgebra, Diagonal
using LogDensityProblemsAD: LogDensityProblemsAD
using MacroTools: MacroTools
using ModelingToolkit: ModelingToolkit, parameters, solve
using Optim: Optim
using Parameters: Parameters, @unpack
using Random: Random
using RecursiveArrayTools: RecursiveArrayTools
using Requires: Requires
using SciMLStructures: SciMLStructures
using StanSample: StanSample, SampleModel, read_samples, stan_sample
using TransformVariables: TransformVariables, as, asℝ₊
using TransformedLogDensities: TransformedLogDensities, TransformedLogDensity
using Turing: Turing, InverseGamma, MCMCSerial, MvNormal, NamedDist, Normal,
    logpdf, sample
STANDARD_PROB_GENERATOR(prob, p) = remake(prob; u0 = eltype(p).(prob.u0), p = p)
function STANDARD_PROB_GENERATOR(prob::EnsembleProblem, p)
    return EnsembleProblem(remake(prob.prob; u0 = eltype(p).(prob.prob.u0), p = p))
end

include("turing_inference.jl")
# include("abc_inference.jl")
include("stan_string.jl")
include("stan_inference.jl")
include("dynamichmc_inference.jl")

export turing_inference, stan_inference, dynamichmc_inference
end # module
