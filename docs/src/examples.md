# Bayesian Inference Examples

### Stan

Like in the previous examples, we set up the Lotka-Volterra system and generate
data.

```julia
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

```julia
priors = [truncated(Normal(1.5,0.1),0,2),truncated(Normal(1.0,0.1),0,1.5),
          truncated(Normal(3.0,0.1),0,4),truncated(Normal(1.0,0.1),0,2)]
```

We then give these to the inference function.

```julia
bayesian_result = stan_inference(prob1,t,data,priors;
                                 num_samples=100,num_warmup=500,
                                 vars = (StanODEData(),InverseGamma(4,1)))
```

`InverseGamma(4,1)` is our starting estimation for the variance hyperparameter
of the default `Normal` distribution. The result is a
[Mamba.jl](http://mambajl.readthedocs.io/en/dev/intro.html) chain object.
We can pull out the parameter values via:

```julia
theta1 = bayesian_result.chain_results[:,["theta.1"],:]
theta2 = bayesian_result.chain_results[:,["theta.2"],:]
theta3 = bayesian_result.chain_results[:,["theta.3"],:]
theta4 = bayesian_result.chain_results[:,["theta.4"],:]
```

From these chains we can get our estimate for the parameters via:

```julia
mean(theta1.value[:,:,1])
```

We can get more of a description via:

```julia
Mamba.describe(bayesian_result.chain_results)

# Result

Iterations = 1:100
Thinning interval = 1
Chains = 1,2,3,4
Samples per chain = 100

Empirical Posterior Estimates:
                  Mean         SD        Naive SE        MCSE         ESS
         lp__ -6.15472697 1.657551334 0.08287756670 0.18425029767  80.9314979
accept_stat__  0.90165904 0.125913744 0.00629568721 0.02781181930  20.4968668
   stepsize__  0.68014975 0.112183047 0.00560915237 0.06468790087   3.0075188
  treedepth__  2.68750000 0.524911975 0.02624559875 0.10711170182  24.0159141
 n_leapfrog__  6.77000000 4.121841086 0.20609205428 0.18645821695 100.0000000
  divergent__  0.00000000 0.000000000 0.00000000000 0.00000000000         NaN
     energy__  9.12245750 2.518330231 0.12591651153 0.32894488320  58.6109941
     sigma1.1  0.57164997 0.128579363 0.00642896816 0.00444242658 100.0000000
     sigma1.2  0.58981422 0.131346442 0.00656732209 0.00397310122 100.0000000
       theta1  1.50237077 0.008234095 0.00041170473 0.00025803930 100.0000000
       theta2  0.99778276 0.009752574 0.00048762870 0.00009717115 100.0000000
       theta3  3.00087782 0.009619775 0.00048098873 0.00020301023 100.0000000
       theta4  0.99803569 0.008893244 0.00044466218 0.00040886528 100.0000000
      theta.1  1.50237077 0.008234095 0.00041170473 0.00025803930 100.0000000
      theta.2  0.99778276 0.009752574 0.00048762870 0.00009717115 100.0000000
      theta.3  3.00087782 0.009619775 0.00048098873 0.00020301023 100.0000000
      theta.4  0.99803569 0.008893244 0.00044466218 0.00040886528 100.0000000

Quantiles:
                  2.5%        25.0%      50.0%      75.0%       97.5%
         lp__ -10.11994750 -7.0569000 -5.8086150 -4.96936500 -3.81514375
accept_stat__   0.54808912  0.8624483  0.9472840  0.98695850  1.00000000
   stepsize__   0.57975100  0.5813920  0.6440120  0.74276975  0.85282400
  treedepth__   2.00000000  2.0000000  3.0000000  3.00000000  3.00000000
 n_leapfrog__   3.00000000  7.0000000  7.0000000  7.00000000 15.00000000
  divergent__   0.00000000  0.0000000  0.0000000  0.00000000  0.00000000
     energy__   5.54070300  7.2602200  8.7707000 10.74517500 14.91849500
     sigma1.1   0.38135240  0.4740865  0.5533195  0.64092575  0.89713635
     sigma1.2   0.39674703  0.4982615  0.5613655  0.66973025  0.88361407
       theta1   1.48728600  1.4967650  1.5022750  1.50805500  1.51931475
       theta2   0.97685115  0.9914630  0.9971435  1.00394250  1.01765575
       theta3   2.98354100  2.9937575  3.0001450  3.00819000  3.02065950
       theta4   0.97934128  0.9918495  0.9977415  1.00430750  1.01442975
      theta.1   1.48728600  1.4967650  1.5022750  1.50805500  1.51931475
      theta.2   0.97685115  0.9914630  0.9971435  1.00394250  1.01765575
      theta.3   2.98354100  2.9937575  3.0001450  3.00819000  3.02065950
      theta.4   0.97934128  0.9918495  0.9977415  1.00430750  1.01442975
```

More extensive information about the distributions is given by the plots:

```julia
plot_chain(bayesian_result)
```

### Turing

This case we will build off of the Stan example. Note that `turing_inference`
does not require the use of the `@ode_def` macro like Stan does, but it will
still work with macro-defined functions. Thus, using the same setup as before,
we simply give the setup to:

```julia
bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500)
```

The result is a [`MCMCChains.jl`](https://github.com/TuringLang/MCMCChains.jl)
chains object. The chain for the first parameter is then given by:

```julia
bayesian_result["theta[1]"]
```

Summary statistics can be also be accessed:
```julia
using StatsBase
describe(bayesian_result)
```

The chain can be analysed by the trace plots and other plots obtained by:

```julia
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

```julia
bayesian_result_hmc = dynamichmc_inference(prob1, Tsit5(), t, data, [Normal(1.5, 1)], [bridge(ℝ, ℝ⁺, )])
```

A tuple with summary statistics and the chain values is returned.
The chain for the `i`th parameter is given by:

```julia
bayesian_result_hmc[1][i]
```
For accessing the various summary statistics:

```julia
DynamicHMC.NUTS_statistics(bayesian_result_dynamic[2])
```
Some details about the NUTS sampler can be obtained from:

```julia
bayesian_result_dynamic[3]
```

In case of `dynamic_inference` the trace plots for the `i`th parameter can be obtained by:

```julia
plot(bayesian_result_hmc[1][i])
```

For a better idea of the summary statistics and plotting you can take a look at the [benchmarks](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl).
