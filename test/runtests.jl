using SafeTestsets
const LONGER_TESTS = false

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV,"APPVEYOR")
const is_TRAVIS = haskey(ENV,"TRAVIS")

@time begin
if GROUP == "All" || GROUP == "Core"
    @time @testset "DynamicHMC" begin include("dynamicHMC.jl") end
    @time @testset "Stan_String" begin include("stan_string.jl") end
    @time @testset "Turing" begin include("turing.jl") end # Doesn't work on v0.6
    @time @testset "ABC" begin include("abc.jl") end
end
end

if GROUP == "Stan"
    if is_TRAVIS
      using Pkg
      Pkg.add("CmdStan")
    end
    @time @testset "Stan" begin include("stan.jl") end
end
