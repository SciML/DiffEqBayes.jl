# Bayesian Inference of ODE

For this tutorial, we will show how to do Bayesian inference to infer the parameters of
the Lotka-Volterra equations using each of the three backends:

  - Turing.jl
  - Stan.jl
  - DynamicHMC.jl

## Setup

First, let's set up our ODE and the data. For the data, we will simply solve the ODE and
take that solution at some known parameters as the dataset. This looks like the following:

```@example all
using DiffEqBayes, ParameterizedFunctions, OrdinaryDiffEq, RecursiveArrayTools,
      Distributions
f1 = @ode_def LotkaVolterra begin
    dx = a * x - x * y
    dy = -3 * y + x * y
end a

p = [1.5]
u0 = [1.0, 1.0]
tspan = (0.0, 10.0)
prob1 = ODEProblem(f1, u0, tspan, p)

σ = 0.01                         # noise, fixed for now
t = collect(1.0:10.0)   # observation times
sol = solve(prob1, Tsit5())
priors = [Normal(1.5, 1)]
randomized = VectorOfArray([(sol(t[i]) + σ * randn(2)) for i in 1:length(t)])
data = convert(Array, randomized)
```

## Inference Methods

### Stan

```@example all
using CmdStan #required for using the Stan backend
bayesian_result_stan = stan_inference(prob1, t, data, priors)
```

### Turing

```@example all
bayesian_result_turing = turing_inference(prob1, Tsit5(), t, data, priors)
```

### DynamicHMC

We can use [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) as the backend
for sampling with the `dynamic_inference` function. It is similarly used as follows:

```@example all
bayesian_result_hmc = dynamichmc_inference(prob1, Tsit5(), t, data, priors)
```

## More Information

For a better idea of the summary statistics and plotting, you can take a look at the [benchmarks](https://github.com/SciML/SciMLBenchmarks.jl).
