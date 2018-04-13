using DiffEqBayes
using Base.Test

tic()
@time @testset "DynamicHMC" begin include("dynamicHMC.jl") end
@time @testset "Stan_String" begin include("stan_string.jl") end
@time @testset "Stan" begin include("stan.jl") end
@time @testset "Turing" begin include("turing.jl") end # Doesn't work on v0.6
@time @testset "ABC" begin include("abc.jl") end
toc()
