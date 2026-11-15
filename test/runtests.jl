using SciMLTesting

if get(ENV, "GROUP", "All") == "Stan"
    using SafeTestsets: @safetestset
    stan_dir = joinpath(@__DIR__, "Stan")
    for file in sort!(filter(f -> endswith(f, ".jl"), readdir(stan_dir)))
        @eval @safetestset $file begin
            include(joinpath(@__DIR__, "Stan", $file))
        end
    end
else
    run_tests()
end
