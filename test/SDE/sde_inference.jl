using DiffEqBayes
using DiffEqBase: SDEProblem
using Distributions
using Test

@testset "sde_inference basic" begin
    f!(du, u, p, t) = (du[1] = p[1] * u[1])
    g!(du, u, p, t) = (du[1] = 0.1)
    prob = SDEProblem(f!, g!, [1.0], (0.0, 1.0), [0.5])
    data = rand(1, 10)
    t = collect(range(0.0, 1.0, length = 10))
    priors = [Normal(0.5, 0.1)]
    chain = sde_inference(prob, data, t, priors; num_samples = 10)
    @test chain !== nothing
end
