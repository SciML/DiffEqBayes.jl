using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions, RecursiveArrayTools
using Parameters, Distributions, BlackBoxOptim, Base.Test

pf_func = function (du,u,p,t)
    du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
    du[2] = -3.0 * u[2] + u[1]*u[2]
  end
  u0 = [1.0;1.0]
  tspan = (0.0,10.0)
  p = [1.5,1.0]
  prob1 = ODEProblem(pf_func,u0,tspan,p)
  sol = solve(prob1,Tsit5())
  t = collect(linspace(0,10,200))
  function generate_data(sol,t)
    randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
    data = convert(Array,randomized)
  end
  aggregate_data = convert(Array,VectorOfArray([generate_data(sol,t) for i in 1:100]))
  
  distributions = [fit_mle(Normal,aggregate_data[i,j,:]) for i in 1:2, j in 1:200]
  priors = [Truncated(Normal(1.5,0.1),0,2),Truncated(Normal(1.0,0.1),0,1.5)]
  obj = map_inference(prob1,Tsit5(),t,distributions,priors;maxiters=10000,verbose=false)
  bound1 = Tuple{Float64, Float64}[(0.5, 5),(0.5, 5)]
  result = bboptimize(p->obj(p);SearchRange = bound1, MaxSteps = 11e3)
  @test result.archive_output.best_candidate ≈ [1.5,1.0] atol = 1e-1

  