# Bayesian Inference Examples

### Stan

Like in the previous examples, we set up the Lotka-Volterra system and generate
data.

```@example all
using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions,
      RecursiveArrayTools, Distributions, Test

f1 = @ode_def begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
p = [1.5,1.0,3.0,1.0]
u0 = [1.0,1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,p)
sol = solve(prob1,Tsit5())
t = collect(range(1,stop=10,length=10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
```

Here we now give Stan an array of prior distributions for our parameters. Since
the parameters of our differential equation must be positive, we utilize
truncated Normal distributions to make sure that is satisfied in the result:

```@example all
priors = [truncated(Normal(1.5,0.1),0,2),truncated(Normal(1.0,0.1),0,1.5),
          truncated(Normal(3.0,0.1),0,4),truncated(Normal(1.0,0.1),0,2)]
```

We then give these to the inference function.

```@example all
bayesian_result = stan_inference(prob1,t,data,priors;
                                 num_samples=100,num_warmup=500,
                                 vars = (StanODEData(),InverseGamma(4,1)))
```

`InverseGamma(4,1)` is our starting estimation for the variance hyperparameter
of the default `Normal` distribution. The result is a
[Mamba.jl](http://mambajl.readthedocs.io/en/dev/intro.html) chain object.
We can pull out the parameter values via:

```@example all
theta1 = bayesian_result.chain_results[:,["theta.1"],:]
theta2 = bayesian_result.chain_results[:,["theta.2"],:]
theta3 = bayesian_result.chain_results[:,["theta.3"],:]
theta4 = bayesian_result.chain_results[:,["theta.4"],:]
```

From these chains we can get our estimate for the parameters via:

```@example all
mean(theta1.value[:,:,1])
```

We can get more of a description via:

```@example all
Mamba.describe(bayesian_result.chain_results)
```

More extensive information about the distributions is given by the plots:

```@example all
plot_chain(bayesian_result)
```

### Turing

This case we will build off of the Stan example. Note that `turing_inference`
does not require the use of the `@ode_def` macro like Stan does, but it will
still work with macro-defined functions. Thus, using the same setup as before,
we simply give the setup to:

```@example all
bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500)
```

The result is a [`MCMCChains.jl`](https://github.com/TuringLang/MCMCChains.jl)
chains object. The chain for the first parameter is then given by:

```@example all
bayesian_result["theta[1]"]
```

Summary statistics can be also be accessed:

```@example all
using StatsBase
describe(bayesian_result)
```

The chain can be analysed by the trace plots and other plots obtained by:

```@example all
using StatsPlots
plot(bayesian_result)
```

### DynamicHMC

We can use [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) as the backend
for sampling with the `dynamic_inference` function. It supports any `DEProblem`,
`priors` can be passed as an array of [Distributions.jl](https://juliastats.github.io/Distributions.jl/dev/)
distributions, passing `initial` values is optional and in case where the user has a firm understanding of the
domain the parameter values will lie in, `transformations` can be used to pass an array of constraints for the parameters
as an array of [Transformations](https://github.com/tpapp/ContinuousTransformations.jl).

```@example all
bayesian_result_hmc = dynamichmc_inference(prob1, Tsit5(), t, data, [Normal(1.5, 1)], [bridge(ℝ, ℝ⁺, )])
```

A tuple with summary statistics and the chain values is returned.
The chain for the `i`th parameter is given by:

```@example all
bayesian_result_hmc[1][i]
```
For accessing the various summary statistics:

```@example all
DynamicHMC.NUTS_statistics(bayesian_result_dynamic[2])
```
Some details about the NUTS sampler can be obtained from:

```@example all
bayesian_result_dynamic[3]
```

In case of `dynamic_inference` the trace plots for the `i`th parameter can be obtained by:

```@example all
plot(bayesian_result_hmc[1][i])
```

For a better idea of the summary statistics and plotting you can take a look at the [benchmarks](https://github.com/SciML/SciMLBenchmarks.jl).
