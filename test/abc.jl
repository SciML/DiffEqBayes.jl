using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions, Distances, StatsBase, RecursiveArrayTools
using Test

# One parameter case
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

bayesian_result = abc_inference(prob1,Tsit5(),t,data,priors;
                                num_samples=500,ϵ = 0.001)

@test mean(bayesian_result.parameters, weights(bayesian_result.weights)) ≈ 1.5 atol=0.1

# custom distance-function
weights_ = ones(size(data)) # weighted data
for i = 1:3:length(data)
    weights_[i] = 0
    data[i] = 1e20 # to test that those points are indeed not used
end
distfn = function (d1, d2)
    d = 0.0
    for i=1:length(d1)
        d += (d1[i] - d2[i])^2 * weights_[i]
    end
    return sqrt(d)
end
bayesian_result = abc_inference(prob1,Tsit5(),t,data,priors;
                                num_samples=500, ϵ = 0.001,
                                distancefunction = distfn)

@test mean(bayesian_result.parameters, weights(bayesian_result.weights)) ≈ 1.5 atol=0.1

# Four parameter case
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
priors = [Truncated(Normal(1.5,1),0,2),Truncated(Normal(1.0,1),0,1.5),
          Truncated(Normal(3.0,1),0,4),Truncated(Normal(1.0,1),0,2)]

bayesian_result = abc_inference(prob1,Tsit5(),t,data,priors; num_samples=500)

meanvals = mean(bayesian_result.parameters, weights(bayesian_result.weights), 1)

@test meanvals[1] ≈ 1.5 atol=3e-1
@test meanvals[2] ≈ 1.0 atol=3e-1
@test meanvals[3] ≈ 3.0 atol=3e-1
@test meanvals[4] ≈ 1.0 atol=3e-1
