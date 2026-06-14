using SafeTestsets

@time @safetestset "DynamicHMC" begin
    include("dynamicHMC.jl")
end
@time @safetestset "Turing" begin
    include("turing.jl")
end
# @time @safetestset "ABC" begin include("abc.jl") end
