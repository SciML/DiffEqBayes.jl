@setup_workload begin
    precompile_problem = DiffEqBase.ODEProblem(
        (u, p, t) -> p[1] * u,
        [1.0],
        (0.0, 1.0),
        [0.1]
    )
    precompile_prior = Normal(0.0, 1.0)

    @compile_workload begin
        STANDARD_PROB_GENERATOR(precompile_problem, [0.1])
        stan_string(precompile_prior)
        StanODEData()
        StanResult(:model, 0, (;))
    end
end
