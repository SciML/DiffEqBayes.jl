using SafeTestsets
using SciMLTesting

const GROUP = current_group()

if GROUP == "JET"
    # The jet env is activated and instantiated but the repo root is NOT
    # `Pkg.develop`ed (matching the original branch), so this keeps the explicit
    # dispatch rather than run_tests' env mechanism (which always develops).
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "jet"))
    Pkg.instantiate()
    @time @safetestset "JET" begin
        include("jet/jet_tests.jl")
    end
else
    run_tests(;
        core = joinpath(@__DIR__, "core_tests.jl"),
        groups = Dict(
            "Stan" => joinpath(@__DIR__, "stan_tests.jl"),
        ),
        all = ["Core", "Stan"],
    )
end
