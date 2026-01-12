using SafeTestsets
const LONGER_TESTS = false

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    @time @safetestset "DynamicHMC" begin
        include("dynamicHMC.jl")
    end
    @time @safetestset "Turing" begin
        include("turing.jl")
    end
    # @time @safetestset "ABC" begin include("abc.jl") end
end

if GROUP == "Stan" || GROUP == "All"
    @time @safetestset "Stan_String" begin
        include("stan_string.jl")
    end
    @time @safetestset "Stan" begin
        include("stan.jl")
    end
end

if GROUP == "JET"
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "jet"))
    Pkg.instantiate()
    @time @safetestset "JET" begin
        include("jet/jet_tests.jl")
    end
end
