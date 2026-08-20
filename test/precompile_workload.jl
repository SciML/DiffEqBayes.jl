using DiffEqBayes
using Distributions: Normal
using Test

@testset "Precompile workload API" begin
    @test DiffEqBayes.stan_string(Normal(0.0, 1.0)) == "normal(0.0, 1.0)"
    @test StanODEData() isa StanODEData
    @test StanResult(:model, 0, (;)).return_code == 0
end
