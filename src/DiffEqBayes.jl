module DiffEqBayes
using DiffEqBase, CmdStan, Distributions, Turing, MacroTools, Mamba
using ParameterizedFunctions, RecursiveArrayTools
#using DynamicHMC, TransformVariables, LogDensityProblems
using Parameters, Distributions, Optim, Requires
using Distances, ApproxBayes, StatsPlots

STANDARD_PROB_GENERATOR(prob,p) = remake(prob;u0=eltype(p).(prob.u0),p=p)
STANDARD_PROB_GENERATOR(prob::EnsembleProblem,p) = EnsembleProblem(remake(prob.prob;u0=eltype(p).(prob.prob.u0),p=p))

include("turing_inference.jl")
include("abc_inference.jl")

function __init__()
    @require CmdStan="593b3428-ca2f-500c-ae53-031589ec8ddd" begin
        include("stan_inference.jl")
        include("stan_string.jl")
        export StanModel, stan_inference, stan_string, StanODEData
    end

    @require DynamicHMC="bbc10e6e-7c05-544b-b16e-64fede858acb" begin
        include("dynamichmc_inference.jl")
        export dynamichmc_inference
    end

    @require StatsPlots="f3b207a7-027a-5e70-b257-86293d7955fd" begin
        include("utils.jl")
        export plot_chain
    end
end

export turing_inference, abc_inference

end # module
