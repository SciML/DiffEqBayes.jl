using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions,
      RecursiveArrayTools, Distributions, Test

println("One parameter case")
f1 = @ode_def begin
  dx = a*x - x*y
  dy = -3y + x*y
end a
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5]
prob1 = ODEProblem(f1,u0,tspan,p)
sol = solve(prob1,Tsit5())
t = collect(range(1,stop=10,length=10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = [Truncated(Normal(1.5,0.1),0,2)]

bayesian_result = stan_inference(prob1,t,data,priors;num_samples=300,
                                 num_warmup=500,likelihood=Normal,
                                 vars =(StanODEData(),InverseGamma(3,2)))

theta1 = bayesian_result.chain_results[:,["theta.1"],:]
@test mean(theta1.value) ≈ 1.5 atol=3e-1


println("Four parameter case")
f1 = @ode_def begin
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
priors = [Truncated(Normal(1.5,0.01),0,2),Truncated(Normal(1.0,0.01),0,1.5),
          Truncated(Normal(3.0,0.01),0,4),Truncated(Normal(1.0,0.01),0,2)]

bayesian_result = stan_inference(prob1,t,data,priors;num_samples=100,num_warmup=500,vars =(StanODEData(),InverseGamma(4,1)))
theta1 = bayesian_result.chain_results[:,["theta.1"],:]
theta2 = bayesian_result.chain_results[:,["theta.2"],:]
theta3 = bayesian_result.chain_results[:,["theta.3"],:]
theta4 = bayesian_result.chain_results[:,["theta.4"],:]
@test mean(theta1.value[:,:,1]) ≈ 1.5 atol=1e-1
@test mean(theta2.value[:,:,1]) ≈ 1.0 atol=1e-1
@test mean(theta3.value[:,:,1]) ≈ 3.0 atol=1e-1
@test mean(theta4.value[:,:,1]) ≈ 1.0 atol=1e-1
