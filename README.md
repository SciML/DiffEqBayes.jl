# DiffEqBayes.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/DiffEqBayes.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DiffEqBayes.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/DiffEqBayes.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/DiffEqBayes.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/DiffEqBayes.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/DiffEqBayes.jl?branch=master)

This repository is a set of extension functionality for estimating the parameters of differential equations using Bayesian methods. It allows the choice of using [Stan.jl](https://github.com/goedman/Stan.jl), [Turing.jl](https://github.com/yebai/Turing.jl), [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) and [ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl) to perform a Bayesian estimation of a differential equation problem specified via the [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) interface.

To begin you first need to add this repository using the following command.
```julia
Pkg.add("DiffEqBayes")
using DiffEqBayes
```

### stan_inference

```julia
stan_inference(prob::DiffEqBase.DEProblem,t,data,priors = nothing;
               alg=:rk45, num_samples=1000, num_warmup=1000, 
               reltol=1e-3, abstol=1e-6, maxiter=Int(1e5),likelihood=Normal,
               vars=(StanODEData(),InverseGamma(3,3)),nchains=1, sample_u0 = false, 
               save_idxs = nothing, diffeq_string = nothing, printsummary = true)
```

`stan_inference` uses [CmdStan.jl](https://github.com/StanJulia/CmdStan.jl)
to perform the Bayesian inference. The
[Stan installation process](https://stanjulia.github.io/CmdStan.jl/stable/INSTALLATION.html)
is required to use this function. The first argument is a `DEProblem`, `t` is the array of time
and `data` is the array where the first dimension (columns) corresponds to the
array of system values. `priors` is an array of prior distributions for each
parameter, specified via a [Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/)
type. `alg` is a choice between `:rk45` and `:bdf`, the two internal integrators
of Stan. `num_samples` is the number of samples to take per chain, and `num_warmup`
is the number of MCMC warmup steps. `abstol` and `reltol` are the keyword
arguments for the internal integrator. `likelihood` is the likelihood distribution
to use with the arguments from `vars`, and `vars` is a tuple of priors for the
distributions of the likelihood hyperparameters. The special value `StanODEData()`
in this tuple denotes the position that the ODE solution takes in the likelihood's
parameter list. With the `diffeq_string` kwarg you can pass in a complex ODE specification
if the need arises. 

### turing_inference

```julia
turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors;
                    likelihood_dist_priors = [InverseGamma(2, 3)], 
                    likelihood = (u,p,t,σ) -> MvNormal(u, σ[1]*ones(length(u))),
                    num_samples=1000, sampler = Turing.NUTS(0.65),
                    syms = [Turing.@varname(theta[i]) for i in 1:length(priors)],
                    sample_u0 = false, save_idxs = nothing, progress = false, kwargs...)
```

`turing_inference` uses [Turing.jl](https://github.com/TuringLang/Turing.jl) to
perform its parameter inference. `prob` can be any `DEProblem` with a corresponding
`alg` choice. `t` is the array of time points and `data` is the set of
observations for the differential equation system at time point `t[i]` (or higher
dimensional). `priors` is an array of prior distributions for each
parameter, specified via a
[Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/)
type. `num_samples` is the number of samples per MCMC chain. The extra `kwargs` are given to the internal differential
equation solver.

### dynamichmc_inference

```julia
dynamichmc_inference(problem::DiffEqBase.DEProblem, algorithm, t, data,parameter_priors, 
                    parameter_transformations=as(Vector, asℝ₊, length(parameter_priors));
                    σ_priors = fill(Normal(0, 5), size(data, 1)),
                    rng = Random.GLOBAL_RNG, num_samples = 1000,
                    AD_gradient_kind = Val(:ForwardDiff),solve_kwargs = (), 
                    mcmc_kwargs = (initialization = (q = zeros(length(parameter_priors) + 2),),), sample_u0 = false)
```

`dynamichmc_inference` uses [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) to
 perform the bayesian parameter estimation. `prob` can be any `DEProblem`, `data` is the set
 of observations for our model which is to be used in the Bayesian Inference process. `priors` represent the
 choice of prior distributions for the parameters to be determined, passed as an array of [Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/) 
 distributions. `t` is the array of time points. `parameter_transformations` is an array of [Tranformations](https://github.com/tpapp/ContinuousTransformations.jl) 
 imposed for constraining the parameter values to specific domains. `rng` is the random number generator used for MCMC. Defaults to the global one.
`num_samples` is the number of MCMC draws (default: 1000) `AD_gradient_kind` is passed on to `LogDensityProblems.ADgradient`, 
 make sure to `import`the corresponding library. `solve_kwargs` is passed on to `solve`
`mcmc_kwargs` are passed on as keyword arguments to `DynamicHMC.mcmc_with_warmup`


### abc_inference

```julia
abc_inference(prob::DEProblem, alg, t, data, priors; ϵ=0.001,
     distancefunction = euclidean, ABCalgorithm = ABCSMC, progress = false,
     num_samples = 500, maxiterations = 10^5, kwargs...)
```

`abc_inference` uses [ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl) which uses Approximate Bayesian Computation (ABC) to
perform its parameter inference. `prob` can be any `DEProblem` with a corresponding
`alg` choice. `t` is the array of time points and `data[:,i]` is the set of
observations for the differential equation system at time point `t[i]` (or higher
dimensional). `priors` is an array of prior distributions for each
parameter, specified via a
[Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/)
type. `num_samples` is the number of posterior samples. `ϵ` is the target
distance between the data and simulated data. `distancefunction` is a distance metric specified from the
[Distances.jl](https://github.com/JuliaStats/Distances.jl)
package, the default is `euclidean`. `ABCalgorithm` is the ABC algorithm to use, options are `ABCSMC` or `ABCRejection` from
[ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl), the default
is the former which is more efficient. `maxiterations` is the maximum number of iterations before the algorithm terminates. The extra `kwargs` are given to the internal differential
equation solver.


 ## Example

 ```julia
 using ParameterizedFunctions, OrdinaryDiffEq, RecursiveArrayTools, Distributions
 f1 = @ode_def LotkaVolterra begin
  dx = a*x - x*y
  dy = -3*y + x*y
 end a

 p = [1.5]
 u0 = [1.0,1.0]
 tspan = (0.0,10.0)
 prob1 = ODEProblem(f1,u0,tspan,p)

 σ = 0.01                         # noise, fixed for now
 t = collect(1.:10.)   # observation times
 sol = solve(prob1,Tsit5())
 priors = [Normal(1.5, 1)]
 randomized = VectorOfArray([(sol(t[i]) + σ * randn(2)) for i in 1:length(t)])
 data = convert(Array,randomized)
 
 using CmdStan #required for using the Stan backend
 bayesian_result_stan = stan_inference(prob1,t,data,priors)

 bayesian_result_turing = turing_inference(prob1,Tsit5(),t,data,priors)
 
 using DynamicHMC
 bayesian_result_hmc = dynamichmc_inference(prob1, Tsit5(), t, data, priors)

 bayesian_result_abc = abc_inference(prob1, Tsit5(), t, data, priors)
```
### Using save_idxs to declare observables

```julia
 sol = solve(prob1,Tsit5(),save_idxs=[1])
 randomized = VectorOfArray([(sol(t[i]) + σ * randn(1)) for i in 1:length(t)])
 data = convert(Array,randomized)

 using CmdStan #required for using the Stan backend
 bayesian_result_stan = stan_inference(prob1,t,data,priors,save_idxs=[1])

 bayesian_result_turing = turing_inference(prob1,Tsit5(),t,data,priors,save_idxs=[1])
 
 using DynamicHMC
 bayesian_result_hmc = dynamichmc_inference(prob1,Tsit5(),t,data,priors,solve_kwargs = (save_idxs = [1],))

 bayesian_result_abc = abc_inference(prob1,Tsit5(),t,data,priors,save_idxs=[1])
 ```
