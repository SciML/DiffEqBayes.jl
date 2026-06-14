using SafeTestsets

@time @safetestset "Stan_String" begin
    include("stan_string.jl")
end
@time @safetestset "Stan" begin
    include("stan.jl")
end
