using DiffEqBayes
using Test

@time begin
@time @testset "Turing" begin include("turing.jl") end
@time @testset "DynamicHMC" begin include("dynamicHMC.jl") end
@time @testset "Stan_String" begin include("stan_string.jl") end
@time @testset "Stan" begin include("stan.jl") end
@time @testset "ABC" begin include("abc.jl") end
end
