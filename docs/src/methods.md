# Bayesian Methods

The following methods require DiffEqBayes.jl:

```julia
using Pkg
Pkg.add("DiffEqBayes")
using DiffEqBayes
```

### stan_inference

```julia
stan_inference(prob::ODEProblem, t, data, priors = nothing; alg = :rk45,
               num_samples = 1000, num_warmups = 1000, reltol = 1e-3,
               abstol = 1e-6, maxiter = Int(1e5), likelihood = Normal,
               vars = (StanODEData(), InverseGamma(2, 3)))
```

`stan_inference` uses [Stan.jl](https://stanjulia.github.io/CmdStan.jl/latest/INTRO/)
to perform the Bayesian inference. The
[Stan installation process](https://stanjulia.github.io/CmdStan.jl/latest/INSTALLATION/)
is required to use this function. `t` is the array of time
and `data` is the array where the first dimension (columns) corresponds to the
array of system values. `priors` is an array of prior distributions for each
parameter, specified via a [Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/)
type. `alg` is a choice between `:rk45` and `:bdf`, the two internal integrators
of Stan. `num_samples` is the number of samples to take per chain, and `num_warmups`
is the number of MCMC warm-up steps. `abstol` and `reltol` are the keyword
arguments for the internal integrator. `likelihood` is the likelihood distribution
to use with the arguments from `vars`, and `vars` is a tuple of priors for the
distributions of the likelihood hyperparameters. The special value `StanODEData()`
in this tuple denotes the position that the ODE solution takes in the likelihood's
parameter list.

### turing_inference

```julia
function turing_inference(prob::DiffEqBase.DEProblem, alg, t, data, priors;
                          likelihood_dist_priors, likelihood, num_samples = 1000,
                          sampler = Turing.NUTS(num_samples, 0.65), syms, kwargs...)
end
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
dynamichmc_inference(prob::DEProblem, alg, t, data, priors, transformations;
                     σ = 0.01, ϵ = 0.001, initial = Float64[])
```

`dynamichmc_inference` uses [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) to
perform the Bayesian parameter estimation. `prob` can be any `DEProblem`, `data` is the set
of observations for our model which is to be used in the Bayesian Inference process. `priors` represent the
choice of prior distributions for the parameters to be determined, passed as an array of [Distributions.jl]
(https://juliastats.github.io/Distributions.jl/dev/) distributions. `t` is the array of time points. `transformations`
is an array of [Transformations](https://github.com/tpapp/ContinuousTransformations.jl) imposed for constraining the
parameter values to specific domains. `initial` values for the parameters can be passed, if not passed the means of the
`priors` are used. `ϵ` can be used as a kwarg to pass the initial step size for the NUTS algorithm.

### abc_inference

```julia
abc_inference(prob::DEProblem, alg, t, data, priors; ϵ = 0.001,
              distancefunction = euclidean, ABCalgorithm = ABCSMC, progress = false,
              num_samples = 500, maxiterations = 10^5, kwargs...)
```

`abc_inference` uses [ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl), which uses Approximate Bayesian Computation (ABC) to
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
