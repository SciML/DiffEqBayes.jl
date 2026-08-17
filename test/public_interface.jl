using DiffEqBayes
using Distributions: Normal
using OrdinaryDiffEq: Tsit5
using SciMLBase: ODEProblem, solve
using Test
using Turing

@testset "public Bayesian inference interface" begin
    problem = ODEProblem((du, u, p, t) -> (du[1] = p[1]), [0.0], (0.0, 1.0), [0.0])
    times = [0.5, 1.0]
    observations = zeros(1, length(times))
    solution = solve(problem, Tsit5(); saveat = times)
    @test size(solution, 2) == length(times)
    chain = turing_inference(
        problem, Tsit5(), times, observations, [Normal(0, 1)];
        syms = [:rate], sample_args = (sampler = Turing.MH(), num_samples = 2)
    )

    @test size(chain, 1) > 0
    @test StanODEData() isa StanODEData
    result = StanResult(:model, 0, (;))
    @test result.model === :model
    @test result.return_code === 0
    @test result.chains === (;)
end
