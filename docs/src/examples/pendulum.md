# Bayesian Inference of Pendulum Parameters

In this tutorial, we will perform Bayesian parameter inference of the parameters of a
pendulum.

## Set up simple pendulum problem

```@example pendulum
using DiffEqBayes, OrdinaryDiffEq, RecursiveArrayTools, Distributions, Plots, StatsPlots,
      BenchmarkTools, TransformVariables, CmdStan, DynamicHMC
```

Let's define our simple pendulum problem. Here, our pendulum has a drag term `ω`
and a length `L`.

![pendulum](https://user-images.githubusercontent.com/1814174/59942945-059c1680-942f-11e9-991c-2025e6e4ccd3.jpg)

We get first order equations by defining the first term as the velocity and the
second term as the position, getting:

```@example pendulum
function pendulum(du, u, p, t)
    ω, L = p
    x, y = u
    du[1] = y
    du[2] = -ω * y - (9.8 / L) * sin(x)
end

u0 = [1.0, 0.1]
tspan = (0.0, 10.0)
prob1 = ODEProblem(pendulum, u0, tspan, [1.0, 2.5])
```

## Solve the model and plot

To understand the model and generate data, let's solve and visualize the solution
with the known parameters:

```@example pendulum
sol = solve(prob1, Tsit5())
plot(sol)
```

It's the pendulum, so you know what it looks like. It's periodic, but since we
have not made a small angle assumption, it's not exactly `sin` or `cos`. Because
the true dampening parameter `ω` is 1, the solution does not decay over time,
nor does it increase. The length `L` determines the period.

## Create some dummy data to use for estimation

We now generate some dummy data to use for estimation

```@example pendulum
t = collect(range(1, stop = 10, length = 10))
randomized = VectorOfArray([(sol(t[i]) + 0.01randn(2)) for i in 1:length(t)])
data = convert(Array, randomized)
```

Let's see what our data looks like on top of the real solution

```@example pendulum
scatter!(data')
```

This data captures the non-dampening effect and the true period, making it
perfect for attempting a Bayesian inference.

## Perform Bayesian Estimation

Now let's fit the pendulum to the data. Since we know our model is correct,
this should give us back the parameters that we used to generate the data!
Define priors on our parameters. In this case, let's assume we don't have much
information, but have a prior belief that ω is between 0.1 and 3.0, while the
length of the pendulum L is probably around 3.0:

```@example pendulum
priors = [
    truncated(Normal(0.1, 1.0), lower = 0.0),
    truncated(Normal(3.0, 1.0), lower = 0.0)
]
```

Finally, let's run the estimation routine from DiffEqBayes.jl with the Turing.jl backend to check if we indeed recover the parameters!

```@example pendulum
bayesian_result = turing_inference(prob1, Tsit5(), t, data, priors;
    syms = [:omega, :L], sample_args = (num_samples = 10_000,))
```

Notice that while our guesses had the wrong means, the learned parameters converged
to the correct means, meaning that it learned good posterior distributions for the
parameters. To look at these posterior distributions on the parameters, we can
examine the chains:

```@example pendulum
plot(bayesian_result)
```

As a diagnostic, we will also check the parameter chains. The chain is the MCMC
sampling process. The chain should explore parameter space and converge reasonably
well, and we should be taking a lot of samples after it converges (it is these
samples that form the posterior distribution!)

```@example pendulum
plot(bayesian_result, colordim = :parameter)
```

Notice that after a while these chains converge to a “fuzzy line”, meaning it
found the area with the most likelihood and then starts to sample around there,
which builds a posterior distribution around the true mean.

DiffEqBayes.jl allows the choice of using Stan.jl, Turing.jl and DynamicHMC.jl
for MCMC, you can also use ApproxBayes.jl for Approximate Bayesian computation algorithms.
Let's compare the timings across the different MCMC backends.
We'll stick with the default arguments and 10,000 samples in each.
However, there is a lot of room for micro-optimization
specific to each package and algorithm combinations,
you might want to do your own experiments for specific problems
to get better understanding of the performance.

```@example pendulum
@btime bayesian_result = turing_inference(prob1, Tsit5(), t, data, priors;
    syms = [:omega, :L], sample_args = (num_samples = 10_000,))
```

```@example pendulum
@btime bayesian_result = stan_inference(prob1, :rk45, t, data, priors;
sample_kwargs = Dict(:num_samples => 10_000), print_summary = false)
```

```@example pendulum
@btime bayesian_result = dynamichmc_inference(prob1, Tsit5(), t, data, priors;
    num_samples = 10_000)
```
