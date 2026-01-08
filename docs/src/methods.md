# Bayesian Methods

The following methods require DiffEqBayes.jl:

```julia
using Pkg
Pkg.add("DiffEqBayes")
using DiffEqBayes
```

### stan_inference

```julia
stan_inference(prob::DiffEqBase.DEProblem, alg, t, data, priors = nothing;
    stanmodel = nothing, likelihood = Normal, vars = (StanODEData(), InverseGamma(3, 3)), sample_u0 = false,
    solve_kwargs = Dict(), diffeq_string = nothing, sample_kwargs = Dict(),
    output_format = :mcmcchains, print_summary = true, tmpdir = mktempdir())
```

`stan_inference` uses [StanSample.jl](https://stanjulia.github.io/StanSample.jl/stable/)
to perform the Bayesian inference. The
[Stan installation process](https://stanjulia.github.io/StanSample.jl/stable/INSTALLATION/)
is required to use this function. Currently `CmdStan v2.34.1` is supported.

`prob` can be any `DEProblem` with a corresponding `alg` choice. `alg` is a choice between `:rk45` and `:bdf`, the two internal integrators of Stan. `t` is the array of time and `data` is the array where the first dimension (columns) corresponds to the array of system values. `priors` is an array of prior distributions for each parameter, specified via a [Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/) type. `likelihood` is the likelihood distribution to use with the arguments from `vars`, and `vars` is a tuple of priors for the distributions of the likelihood hyperparameters. The special value `StanODEData()` in this tuple denotes the position that the ODE solution takes in the likelihood's parameter list.

`solve_kwargs` is a `Dict` and passed to the stan differential equation solver. `solve_kwargs` may contain `save_idxs`, `reltol`, `abstol`, and `maxiter`. `save_idxs` is documented at [`DifferentialEquations.jl`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/). `sample_kwargs` are passed to the stan sampler and accepts `num_samples`, `num_warmups`, `num_cpp_chains` , `num_chains`, `num_threads`, `delta`. Please refer to the [stan documentation for more information](https://mc-stan.org/docs/cmdstan-guide/mcmc-intro.html).

### turing_inference

```julia
turing_inference(prob::DiffEqBase.DEProblem, alg, t, data, priors;
    likelihood_dist_priors, likelihood, syms, sample_u0 = false, progress = false,
    solve_kwargs = Dict(), sample_args = NamedTuple(), sample_kwargs = Dict())
```

`turing_inference` uses [Turing.jl](https://github.com/TuringLang/Turing.jl) to
perform its parameter inference. `prob` can be any `DEProblem` with a corresponding
`alg` choice. `t` is the array of time points and `data` is the set of
observations for the differential equation system at time point `t[i]` (or higher
dimensional). `priors` is an array of prior distributions for each
parameter, specified via a
[Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/)
type.

The `turing_inference` interacts with `SciML.CommonSolve.solve` and `StatsBase.sample`. Both accept many arguments depending on the solver and sampling algorithm.
These arguments are supplied to `turing_inferene` function via `solve_kwargs`, `sample_args`, and `sample_kwargs` arguments. Please refer to [the `solve` documentation](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/) for `solve_kwargs`, e.g. `solve_kwargs = Dict(:save_idxs => [1])`.
The `solve` keyword arguments default to `save_idxs = nothing`. Similarly please refer to [the `sample` documentation](https://turinglang.org/v0.26/docs/using-turing/guide#sampling-multiple-chains) for `sample_args` and `sample_kwargs`. The four positional argument are as following: `sampler`, the sampling algorithm. Sampling from multiple chains is possible serially or parallelly using `parallel_type`. Third `num_samples`, the number of samples per MCMC chain and `n_chains`, the number of MCMC chains. The positional arguments default to the following values.

```julia
sampler = Turing.NUTS(0.65)
parallel_type = MCMCSerial()
num_samples = 1000
n_chains = 1
```

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
