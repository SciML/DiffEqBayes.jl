struct DynamicHMCPosterior
    alg
    problem
    likelihood
    priors
    kwargs
end

function (P::DynamicHMCPosterior)(a)
    @unpack alg, problem, likelihood, priors, kwargs = P

    prob = problem_new_parameters(problem, a)
    sol = solve(prob, alg; kwargs...)
    if any((s.retcode != :Success for s in sol))
        ℓ = -Inf
    else
        ℓ = likelihood(sol)
    end

    if !isfinite(ℓ) && (ℓ ≠ -Inf)
        ℓ = -Inf                # protect against NaN etc, is it needed?
    end
    logpdf_sum = 0
    for i in length(a)
        logpdf_sum += logpdf(priors[i], a[i])
    end
    ℓ + logpdf_sum
end


function dynamichmc_inference(prob::DiffEqBase.DEProblem, alg, t, data, priors, transformations;
                              σ=0.01, ϵ=0.001, initial=Float64[], num_samples=1000,
                              kwargs...)
    likelihood = sol -> sum( sum(logpdf.(Normal(0.0, σ), sol(t) .- data[:, i]))
                             for (i, t) in enumerate(t) )

    dynamichmc_inference(prob, alg, likelihood, priors, transformations;
                         ϵ=ϵ, initial=initial, num_samples=num_samples,
                         kwargs...)
end

function dynamichmc_inference(prob::DiffEqBase.DEProblem, alg, likelihood, priors, transformations;
                              ϵ=0.001, initial=Float64[], num_samples=1000,
                              kwargs...)
    P = DynamicHMCPosterior(alg, prob, likelihood, priors, kwargs)

    transformations_tuple = Tuple(transformations)
    parameter_transformation = TransformationTuple(transformations_tuple) # assuming a > 0
    PT = TransformLogLikelihood(P, parameter_transformation)
    PTG = ForwardGradientWrapper(PT, zeros(length(priors)));

    lower_bound = Float64[]
    upper_bound = Float64[]

    for i in priors
        push!(lower_bound, minimum(i))
        push!(upper_bound, maximum(i))
    end

    # If no initial position is given use local minimum near expectation of priors.
    if length(initial) == 0
        for i in priors
            push!(initial, mean(i))
        end
        initial_opt = Optim.minimizer(optimize(a -> -P(a),lower_bound,upper_bound,initial,Fminbox(GradientDescent())))
    end

    initial_inverse_transformed = Float64[]
    for i in 1:length(initial_opt)
       para = TransformationTuple(transformations[i])
       push!(initial_inverse_transformed,inverse(para, (initial_opt[i], ))[1])
    end
    #println(initial_inverse_transformed)
    sample, NUTS_tuned = NUTS_init_tune_mcmc(PTG,
                                             initial_inverse_transformed,
                                             num_samples, ϵ=ϵ)

    posterior = ungrouping_map(Vector, get_transformation(PT) ∘ get_position, sample)

    return posterior, sample, NUTS_tuned
end
