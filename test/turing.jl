using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions, RecursiveArrayTools
using Test, Distributions
println("One parameter case")
f1 = @ode_def begin
  dx = a*x - x*y
  dy = -3y + x*y
end a
u0 = [1.0,1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,[1.5])
sol = solve(prob1,Tsit5())
t = collect(range(1,stop=10,length=10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = [Normal(1.5,0.01)]

bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500)

@show mean(bayesian_result[:_theta][50:end])

@test mean(bayesian_result[:_theta][50:end])[1] ≈ 1.5 atol=0.1

println("Four parameter case")
f2 = @ode_def begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob2 = ODEProblem(f2,u0,tspan,p)
sol = solve(prob2,Tsit5())
t = collect(range(1,stop=10,length=10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = [Truncated(Normal(1.5,0.01),0,2),Truncated(Normal(1.0,0.01),0,1.5),
          Truncated(Normal(3.0,0.01),0,4),Truncated(Normal(1.0,0.01),0,2)]

bayesian_result = turing_inference(prob2,Tsit5(),t,data,priors;num_samples=500)

@show mean(bayesian_result[:_theta][:,:][50:end])

@test mean(bayesian_result[:_theta][:,1][50:end]) ≈ 1.5 atol=3e-1
@test mean(bayesian_result[:_theta][:,2][50:end]) ≈ 1.0 atol=3e-1
@test mean(bayesian_result[:_theta][:,3][50:end]) ≈ 3.0 atol=3e-1
@test mean(bayesian_result[:_theta][:,4][50:end]) ≈ 1.0 atol=3e-1
