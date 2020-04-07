using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions, RecursiveArrayTools
using DynamicHMC, TransformVariables, LinearAlgebra
using Parameters, Distributions, Optim
using Test

reporter = LogProgressReport(nothing, 100, 60.0) 

f1 = @ode_def LotkaVolterraTest1 begin
  dx = a*x - x*y
  dy = -3*y + x*y
end a
p = [1.5]
u0 = [1.0,1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,p)

σ = 0.01                         # noise, fixed for now
t = collect(range(1,stop=10,length=10))   # observation times
sol = solve(prob1,Tsit5())
randomized = VectorOfArray([(sol(t[i]) + σ * randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
mcmc_kwargs = (initialization = (q = zeros(1 + 2),), reporter = reporter)
bayesian_result = dynamichmc_inference(prob1, Tsit5(), t, data, (Normal(1.5, 1), ),
                                           as(Vector, asℝ₊, 1),mcmc_kwargs=mcmc_kwargs)
@test mean(p.parameters[1] for p in bayesian_result.posterior) ≈ p[1] atol = 0.1

priors = [Normal(1.,0.01),Normal(1.,0.01),Normal(1.5,0.01)]
mcmc_kwargs = (initialization = (q = zeros(3 + 2),), reporter = reporter)
bayesian_result = dynamichmc_inference(prob1, Tsit5(), t, data, priors,
                                           as(Vector, asℝ₊, 3),mcmc_kwargs=mcmc_kwargs,sample_u0=true)

@test mean(p.parameters[1] for p in bayesian_result.posterior) ≈ 1. atol = 0.1
@test mean(p.parameters[2] for p in bayesian_result.posterior) ≈ 1. atol = 0.1
@test mean(p.parameters[3] for p in bayesian_result.posterior) ≈ 1.5 atol = 0.1

sol = solve(prob1,Tsit5(),save_idxs=[1])
randomized = VectorOfArray([(sol(t[i]) + σ * randn(1)) for i in 1:length(t)])
data = convert(Array,randomized)
mcmc_kwargs = (initialization = (q = zeros(1 + 1),), reporter = reporter)

bayesian_result = dynamichmc_inference(prob1, Tsit5(), t, data, (Normal(1.5, 1), ),
                                           as(Vector, asℝ₊, 1),mcmc_kwargs=mcmc_kwargs, save_idxs = [1])
@test mean(p.parameters[1] for p in bayesian_result.posterior) ≈ p[1] atol = 0.1

priors = [Normal(1.,0.001),Normal(1.5,0.001)]
mcmc_kwargs = (initialization = (q = zeros(2 + 1),), reporter = reporter)
bayesian_result = dynamichmc_inference(prob1, Tsit5(), t, data, priors,
                                           as(Vector, asℝ₊, 2),mcmc_kwargs=mcmc_kwargs,save_idxs=[1],sample_u0=true)

@test mean(p.parameters[1] for p in bayesian_result.posterior) ≈ 1. atol = 0.1
@test mean(p.parameters[2] for p in bayesian_result.posterior) ≈ 1.5 atol = 0.1

# With hand-code likelihood function
weights_ = ones(size(data)) # weighted data
for i = 1:3:length(data)
    weights_[i] = 0
    data[i] = 1e20 # to test that those points are indeed not used
end
likelihood = function (sol)
    l = zero(eltype(first(sol)))
    for (i, t) in enumerate(t)
        l += sum(logpdf.(Normal(0.0, σ), sol(t) - data[:, i]) .* weights_[:,i])
    end
    return l
end
@test_broken bayesian_result = dynamichmc_inference(prob1, Tsit5(), likelihood, [Truncated(Normal(1.5, 1), 0, 2)], as((a = asℝ₊,)))
@test_broken mean(bayesian_result[1][1]) ≈ 1.5 atol=1e-1


f1 = @ode_def LotkaVolterraTest4 begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob1 = ODEProblem(f1,u0,tspan,p)
sol = solve(prob1,Tsit5())
t = collect(range(1,stop=10,length=10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = (a = Truncated(Normal(1.5,0.01), 0, 2),
              b = Truncated(Normal(1.0,0.01), 0, 1.5),
              c = Truncated(Normal(3.0,0.01), 0, 4),
              d = Truncated(Normal(1.0, 0.01), 0, 2))
mcmc_kwargs = (initialization = (q = zeros(4 + 2),), reporter = reporter)
bayesian_result = dynamichmc_inference(prob1, Tsit5(), t, data, priors,
          as(Vector, asℝ₊, 4),mcmc_kwargs = mcmc_kwargs)
@test norm(mean([p.parameters for p in bayesian_result.posterior]) .- p, Inf) ≤ 0.1