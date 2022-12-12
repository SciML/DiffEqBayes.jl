var documenterSearchIndex = {"docs":
[{"location":"examples/pendulum/#Bayesian-Inference-of-Pendulum-Parameters","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"","category":"section"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"In this tutorial we will perform Bayesian parameter inference of the parameters of a pendulum.","category":"page"},{"location":"examples/pendulum/#Set-up-simple-pendulum-problem","page":"Bayesian Inference of Pendulum Parameters","title":"Set up simple pendulum problem","text":"","category":"section"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"using DiffEqBayes, OrdinaryDiffEq, RecursiveArrayTools, Distributions, Plots, StatsPlots, BenchmarkTools, TransformVariables, CmdStan, DynamicHMC","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"Let's define our simple pendulum problem. Here our pendulum has a drag term ω and a length L.","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"(Image: pendulum)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"We get first order equations by defining the first term as the velocity and the second term as the position, getting:","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"function pendulum(du,u,p,t)\n    ω,L = p\n    x,y = u\n    du[1] = y\n    du[2] = - ω*y -(9.8/L)*sin(x)\nend\n\nu0 = [1.0,0.1]\ntspan = (0.0,10.0)\nprob1 = ODEProblem(pendulum,u0,tspan,[1.0,2.5])","category":"page"},{"location":"examples/pendulum/#Solve-the-model-and-plot","page":"Bayesian Inference of Pendulum Parameters","title":"Solve the model and plot","text":"","category":"section"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"To understand the model and generate data, let's solve and visualize the solution with the known parameters:","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"sol = solve(prob1,Tsit5())\nplot(sol)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"It's the pendulum, so you know what it looks like. It's periodic, but since we have not made a small angle assumption it's not exactly sin or cos. Because the true dampening parameter ω is 1, the solution does not decay over time, nor does it increase. The length L determines the period.","category":"page"},{"location":"examples/pendulum/#Create-some-dummy-data-to-use-for-estimation","page":"Bayesian Inference of Pendulum Parameters","title":"Create some dummy data to use for estimation","text":"","category":"section"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"We now generate some dummy data to use for estimation","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"t = collect(range(1,stop=10,length=10))\nrandomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])\ndata = convert(Array,randomized)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"Let's see what our data looks like on top of the real solution","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"scatter!(data')","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"This data captures the non-dampening effect and the true period, making it perfect to attempting a Bayesian inference.","category":"page"},{"location":"examples/pendulum/#Perform-Bayesian-Estimation","page":"Bayesian Inference of Pendulum Parameters","title":"Perform Bayesian Estimation","text":"","category":"section"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"Now let's fit the pendulum to the data. Since we know our model is correct, this should give us back the parameters that we used to generate the data! Define priors on our parameters. In this case, let's assume we don't have much information, but have a prior belief that ω is between 0.1 and 3.0, while the length of the pendulum L is probably around 3.0:","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"priors = [Uniform(0.1,3.0), Normal(3.0,1.0)]","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"Finally let's run the estimation routine from DiffEqBayes.jl with the Turing.jl backend to check if we indeed recover the parameters!","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=10_000,\n                                   syms = [:omega,:L])","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"Notice that while our guesses had the wrong means, the learned parameters converged to the correct means, meaning that it learned good posterior distributions for the parameters. To look at these posterior distributions on the parameters, we can examine the chains:","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"plot(bayesian_result)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"As a diagnostic, we will also check the parameter chains. The chain is the MCMC sampling process. The chain should explore parameter space and converge reasonably well, and we should be taking a lot of samples after it converges (it is these samples that form the posterior distribution!)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"plot(bayesian_result, colordim = :parameter)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"Notice that after awhile these chains converge to a \"fuzzy line\", meaning it found the area with the most likelihood and then starts to sample around there, which builds a posterior distribution around the true mean.","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"DiffEqBayes.jl allows the choice of using Stan.jl, Turing.jl and DynamicHMC.jl for MCMC, you can also use ApproxBayes.jl for Approximate Bayesian computation algorithms. Let's compare the timings across the different MCMC backends. We'll stick with the default arguments and 10,000 samples in each since there is a lot of room for micro-optimization specific to each package and algorithm combinations, you might want to do your own experiments for specific problems to get better understanding of the performance.","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"@btime bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;syms = [:omega,:L],num_samples=10_000)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"@btime bayesian_result = stan_inference(prob1,t,data,priors;num_samples=10_000,print_summary=false)","category":"page"},{"location":"examples/pendulum/","page":"Bayesian Inference of Pendulum Parameters","title":"Bayesian Inference of Pendulum Parameters","text":"@btime bayesian_result = dynamichmc_inference(prob1,Tsit5(),t,data,priors;num_samples = 10_000)","category":"page"},{"location":"examples/#Bayesian-Inference-of-ODE","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"","category":"section"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"For this tutorial we will show how to do Bayesian inference to infer the parameters of the Lotka-Volterra equations using each of the three backends:","category":"page"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"Turing.jl\nStan.jl\nDynamicHMC.jl","category":"page"},{"location":"examples/#Setup","page":"Bayesian Inference of ODE","title":"Setup","text":"","category":"section"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"First let's setup our ODE and the data. For the data, we will simply solve the ODE and take that solution at some known parameters as the dataset. This looks like the following:","category":"page"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"using DiffEqBayes, ParameterizedFunctions, OrdinaryDiffEq, RecursiveArrayTools, Distributions\nf1 = @ode_def LotkaVolterra begin\n dx = a*x - x*y\n dy = -3*y + x*y\nend a\n\np = [1.5]\nu0 = [1.0,1.0]\ntspan = (0.0,10.0)\nprob1 = ODEProblem(f1,u0,tspan,p)\n\nσ = 0.01                         # noise, fixed for now\nt = collect(1.:10.)   # observation times\nsol = solve(prob1,Tsit5())\npriors = [Normal(1.5, 1)]\nrandomized = VectorOfArray([(sol(t[i]) + σ * randn(2)) for i in 1:length(t)])\ndata = convert(Array,randomized)","category":"page"},{"location":"examples/#Inference-Methods","page":"Bayesian Inference of ODE","title":"Inference Methods","text":"","category":"section"},{"location":"examples/#Stan","page":"Bayesian Inference of ODE","title":"Stan","text":"","category":"section"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"using CmdStan #required for using the Stan backend\nbayesian_result_stan = stan_inference(prob1,t,data,priors)","category":"page"},{"location":"examples/#Turing","page":"Bayesian Inference of ODE","title":"Turing","text":"","category":"section"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"bayesian_result_turing = turing_inference(prob1,Tsit5(),t,data,priors)","category":"page"},{"location":"examples/#DynamicHMC","page":"Bayesian Inference of ODE","title":"DynamicHMC","text":"","category":"section"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"We can use DynamicHMC.jl as the backend for sampling with the dynamic_inference function. It is similarly used as follows:","category":"page"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"bayesian_result_hmc = dynamichmc_inference(prob1, Tsit5(), t, data, priors)","category":"page"},{"location":"examples/#More-Information","page":"Bayesian Inference of ODE","title":"More Information","text":"","category":"section"},{"location":"examples/","page":"Bayesian Inference of ODE","title":"Bayesian Inference of ODE","text":"For a better idea of the summary statistics and plotting you can take a look at the benchmarks.","category":"page"},{"location":"methods/#Bayesian-Methods","page":"Methods","title":"Bayesian Methods","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"The following methods require DiffEqBayes.jl:","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"]add DiffEqBayes\nusing DiffEqBayes","category":"page"},{"location":"methods/#stan_inference","page":"Methods","title":"stan_inference","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"stan_inference(prob::ODEProblem,t,data,priors = nothing;alg=:rk45,\n               num_samples=1000, num_warmups=1000, reltol=1e-3,\n               abstol=1e-6, maxiter=Int(1e5),likelihood=Normal,\n               vars=(StanODEData(),InverseGamma(2,3)))","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"stan_inference uses Stan.jl to perform the Bayesian inference. The Stan installation process is required to use this function. t is the array of time and data is the array where the first dimension (columns) corresponds to the array of system values. priors is an array of prior distributions for each parameter, specified via a Distributions.jl type. alg is a choice between :rk45 and :bdf, the two internal integrators of Stan. num_samples is the number of samples to take per chain, and num_warmups is the number of MCMC warmup steps. abstol and reltol are the keyword arguments for the internal integrator. likelihood is the likelihood distribution to use with the arguments from vars, and vars is a tuple of priors for the distributions of the likelihood hyperparameters. The special value StanODEData() in this tuple denotes the position that the ODE solution takes in the likelihood's parameter list.","category":"page"},{"location":"methods/#turing_inference","page":"Methods","title":"turing_inference","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"function turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors;\n                              likelihood_dist_priors, likelihood, num_samples=1000,\n                              sampler = Turing.NUTS(num_samples, 0.65), syms, kwargs...)","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"turing_inference uses Turing.jl to perform its parameter inference. prob can be any DEProblem with a corresponding alg choice. t is the array of time points and data is the set of observations for the differential equation system at time point t[i] (or higher dimensional). priors is an array of prior distributions for each parameter, specified via a Distributions.jl type. num_samples is the number of samples per MCMC chain. The extra kwargs are given to the internal differential equation solver.","category":"page"},{"location":"methods/#dynamichmc_inference","page":"Methods","title":"dynamichmc_inference","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"dynamichmc_inference(prob::DEProblem,alg,t,data,priors,transformations;\n                      σ = 0.01,ϵ=0.001,initial=Float64[])","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"dynamichmc_inference uses DynamicHMC.jl to  perform the bayesian parameter estimation. prob can be any DEProblem, data is the set  of observations for our model which is to be used in the Bayesian Inference process. priors represent the  choice of prior distributions for the parameters to be determined, passed as an array of Distributions.jl distributions. t is the array of time points. transformations  is an array of Tranformations imposed for constraining the  parameter values to specific domains. initial values for the parameters can be passed, if not passed the means of the  priors are used. ϵ can be used as a kwarg to pass the initial step size for the NUTS algorithm.","category":"page"},{"location":"methods/#abc_inference","page":"Methods","title":"abc_inference","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"abc_inference(prob::DEProblem, alg, t, data, priors; ϵ=0.001,\n     distancefunction = euclidean, ABCalgorithm = ABCSMC, progress = false,\n     num_samples = 500, maxiterations = 10^5, kwargs...)","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"abc_inference uses ApproxBayes.jl which uses Approximate Bayesian Computation (ABC) to perform its parameter inference. prob can be any DEProblem with a corresponding alg choice. t is the array of time points and data[:,i] is the set of observations for the differential equation system at time point t[i] (or higher dimensional). priors is an array of prior distributions for each parameter, specified via a Distributions.jl type. num_samples is the number of posterior samples. ϵ is the target distance between the data and simulated data. distancefunction is a distance metric specified from the Distances.jl package, the default is euclidean. ABCalgorithm is the ABC algorithm to use, options are ABCSMC or ABCRejection from ApproxBayes.jl, the default is the former which is more efficient. maxiterations is the maximum number of iterations before the algorithm terminates. The extra kwargs are given to the internal differential equation solver.","category":"page"},{"location":"#DiffEqBayes.jl:-Bayesian-Parameter-Estimation-for-Differential-Equations","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"","category":"section"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"This repository is a set of extension functionality for estimating the parameters of differential equations using Bayesian methods. It allows the choice of using CmdStan.jl, Turing.jl, DynamicHMC.jl and ApproxBayes.jl to perform a Bayesian estimation of a differential equation problem specified via the DifferentialEquations.jl interface.","category":"page"},{"location":"#Installation","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"Installation","text":"","category":"section"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"To install DiffEqBayes.jl, use the Julia package manager:","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"using Pkg\nPkg.add(\"DiffEqBayes\")","category":"page"},{"location":"#Contributing","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"Contributing","text":"","category":"section"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"Please refer to the SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages for guidance on PRs, issues, and other matters relating to contributing to SciML.\nSee the SciML Style Guide for common coding practices and other style decisions.\nThere are a few community forums:\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Slack\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Zulip\nOn the Julia Discourse forums\nSee also SciML Community page","category":"page"},{"location":"#Reproducibility","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"</details>","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"</details>","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"using Pkg # hide\nPkg.status(;mode = PKGMODE_MANIFEST) # hide","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"</details>","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"You can also download the \n<a href=\"","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"using TOML\nversion = TOML.parse(read(\"../../Project.toml\",String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\",String))[\"name\"]\nlink = \"https://github.com/SciML/\"*name*\".jl/tree/gh-pages/v\"*version*\"/assets/Manifest.toml\"","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"\">manifest</a> file and the\n<a href=\"","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"using TOML\nversion = TOML.parse(read(\"../../Project.toml\",String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\",String))[\"name\"]\nlink = \"https://github.com/SciML/\"*name*\".jl/tree/gh-pages/v\"*version*\"/assets/Project.toml\"","category":"page"},{"location":"","page":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","title":"DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations","text":"\">project</a> file.","category":"page"}]
}
