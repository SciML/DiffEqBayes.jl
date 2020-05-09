using SafeTestsets
const LONGER_TESTS = false

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV,"APPVEYOR")
const is_TRAVIS = haskey(ENV,"TRAVIS")

@time begin
if GROUP == "All" || GROUP == "Core"
    @time @safetestset "DynamicHMC" begin include("dynamicHMC.jl") end
    @time @safetestset "Turing" begin include("turing.jl") end
    @time @safetestset "ABC" begin include("abc.jl") end
end
end

if GROUP == "Stan"
    if is_TRAVIS
      using Pkg
      Pkg.add("CmdStan")
    end
    @time @safetestset "Stan_String" begin include("stan_string.jl") end
    @time @safetestset "Stan" begin include("stan.jl") end
end
