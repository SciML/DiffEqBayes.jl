using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions, RecursiveArrayTools
using Test, Distributions, SteadyStateDiffEq
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

bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500,
                                   syms=[:a])

@show bayesian_result

@test mean(get(bayesian_result,:a)[1]) ≈ 1.5 atol=3e-1

priors = [Normal(1.,0.01),Normal(1.,0.01),Normal(1.5,0.01)]
bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500,sample_u0 =true,
                                   syms=[:u1,:u2,:a])

@test mean(get(bayesian_result,:a)[1]) ≈ 1.5 atol=3e-1
@test mean(get(bayesian_result,:u1)[1]) ≈ 1.0 atol=3e-1
@test mean(get(bayesian_result,:u2)[1]) ≈ 1.0 atol=3e-1

sol = solve(prob1,Tsit5(),save_idxs=[1])
randomized = VectorOfArray([(sol(t[i]) + .01 * randn(1)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = [Normal(1.5,0.01)]
bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500,
                                   syms=[:a],save_idxs=[1])

@test mean(get(bayesian_result,:a)[1]) ≈ 1.5 atol=3e-1

priors = [Normal(1.,0.01),Normal(1.5,0.01)]
bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500,sample_u0 =true,
                                   syms=[:u1,:a],save_idxs=[1])

@test mean(get(bayesian_result,:a)[1]) ≈ 1.5 atol=3e-1
@test mean(get(bayesian_result,:u1)[1]) ≈ 1.0 atol=3e-1

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
priors = [truncated(Normal(1.5,0.01),0,2),truncated(Normal(1.0,0.01),0,1.5),
          truncated(Normal(3.0,0.01),0,4),truncated(Normal(1.0,0.01),0,2)]

bayesian_result = turing_inference(prob2,Tsit5(),t,data,priors;num_samples=500,
                                   syms = [:a,:b,:c,:d])

@show bayesian_result

@test mean(get(bayesian_result,:a)[1]) ≈ 1.5 atol=3e-1
@test mean(get(bayesian_result,:b)[1]) ≈ 1.0 atol=3e-1
@test mean(get(bayesian_result,:c)[1]) ≈ 3.0 atol=3e-1
@test mean(get(bayesian_result,:d)[1]) ≈ 1.0 atol=3e-1

println("Steady state problem")
function f(du,u,p,t)
    α = p[1]
    du[1] = 2 -  α*u[1]
    du[2] = u[1] - 4u[2]
end

p = [2.0]
u0 = zeros(2)
s_prob = SteadyStateProblem(f,u0,p)
s_sol = solve(s_prob,SSRootfind())
s_sol = solve(s_prob,DynamicSS(Tsit5(),abstol=1e-4,reltol=1e-3))


# true data is 1.00, 0.25
data = [1.05, 0.23]
priors = [truncated(Normal(2.0,0.2),0,3)]
bayesian_result = turing_inference(s_prob,DynamicSS(Tsit5(),abstol=1e-4,reltol=1e-3),
                                   nothing,data,priors;
                                   num_samples=500,
                                   maxiters = 1e6,
                                   syms = [:α])
                                   
@test mean(get(bayesian_result,:α)[1]) ≈ 2.0 atol=3e-1
