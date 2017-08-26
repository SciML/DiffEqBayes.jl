using DiffEqBayes
using Base.Test
using OrdinaryDiffEq, ParameterizedFunctions, RecursiveArrayTools

println("One parameter case")
f1 = @ode_def_nohes LotkaVolterraTest1 begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=1.0 c=3.0 d=1.0
u0 = [1.0,1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan)
sol = solve(prob1,Tsit5())
t = collect(linspace(1,10,10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = [Normal(1.5,1)]

bayesian_result = bayesian_inference(prob1,t,data,priors;num_samples=100,num_warmup=100)
theta1 = bayesian_result.chain_results[:,["theta.1"],:]
@test mean(theta1.value[:,:,1]) ≈ 1.5 atol=1e-1


g1 = @ode_def_nohes LorenzExample begin
  dx = σ*(y-x)
  dy = x*(ρ-z) - y
  dz = x*y - β*z
end σ=10.0 ρ=28.0 β=>2.6666

r0 = [0.1;0.0;0.0]
tspan = (0.0,4.0)
prob = ODEProblem(g1,r0,tspan)
sol = solve(prob,Tsit5())
t = collect(linspace(0.1,4.0,10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(3)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = [Normal(2.6666,1)]

bayesian_result = bayesian_inference(prob,t,data,priors)
theta1 = bayesian_result.chain_results[:,["theta.1"],:]
@test mean(theta1.value[:,:,1]) ≈ 2.6666 atol=1e-1
