using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions, RecursiveArrayTools
using Base.Test
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
priors = [Normal(1.5,0.01)]

bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500)

@show mean(bayesian_result[:theta1][50:end])

@test mean(bayesian_result[:theta1][50:end]) ≈ 1.5 atol=0.1
