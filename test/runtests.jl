using SciMLTesting

# The Stan group is intentionally absent from test_groups.toml: it needs a cmdstan
# build (.github/scripts/stan.sh + JULIA_CMDSTAN_HOME) that the centralized
# grouped-tests.yml matrix cannot provision, so it must not appear as a centralized
# job. Folder-discovery `run_tests()` only accepts groups declared in
# test_groups.toml, so it would reject GROUP="Stan". The dedicated CI.yml workflow
# (which runs stan.sh first) dispatches GROUP="Stan" here under `Pkg.test`, so the
# test environment is already active; include the Stan folder's files directly.
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
