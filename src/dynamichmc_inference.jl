struct DynamicHMCPosterior
    problem
    data
    a_prior
    t
    ϵ_dist
    alg
    kwargs
end

function (P::DynamicHMCPosterior)(θ)
    #println(θ)
    @unpack problem, data, a_prior, t, ϵ_dist, alg, kwargs = P
    a = θ
    try
        prob = problem_new_parameters(problem, a)
        sol = solve(prob, alg; kwargs...)
        ℓ = sum(sum(logpdf.(ϵ_dist, sol(t) .- data[:, i]))
                for (i, t) in enumerate(t))
    catch
        ℓ = -Inf
    end
    if !isfinite(ℓ) && (ℓ ≠ -Inf)
        ℓ = -Inf                # protect against NaN etc, is it needed?
    end
    logpdf_sum = 0
    for i in length(a)
        logpdf_sum += logpdf(a_prior[i], a[i])
    end
    ℓ + logpdf_sum
end

function dynamichmc_inference(prob::DEProblem, alg, t, data, priors, transformations;
                              σ=0.01, ϵ=0.001, initial=Float64[], num_samples=1000, kwargs...)
    P = DynamicHMCPosterior(prob, data, priors, t, Normal(0.0, σ), alg, kwargs)

    transformations_tuple = Tuple(transformations)
    parameter_transformation = TransformationTuple(transformations_tuple) # assuming a > 0
    PT = TransformLogLikelihood(P, parameter_transformation)
    PTG = ForwardGradientWrapper(PT, zeros(length(priors)));

    lower_bound = Float64[]
    upper_bound = Float64[]

    for i in priors
        push!(lower_bound,minimum(i))
        push!(upper_bound,maximum(i))
    end

    if length(initial) == 0
        # If no initial position is given use local minimum near expectation of priors.
        for i in priors
            push!(initial,mean(i))
        end
        initial = Optim.minimizer(optimize(a -> -P(a),initial,lower_bound,upper_bound,Fminbox{GradientDescent}()))
    end

    initial_inverse_transformed = Float64[]
    for i in 1:length(initial)
       para = TransformationTuple(transformations[i])
       push!(initial_inverse_transformed,inverse(para, (initial[i], ))[1])
    end
    #println(initial_inverse_transformed)
    sample, NUTS_tuned = NUTS_init_tune_mcmc(PTG,
                                             initial_inverse_transformed,
                                             num_samples, ϵ=ϵ)

    posterior = ungrouping_map(Vector, get_transformation(PT) ∘ get_position, sample)

    return posterior, sample, NUTS_tuned
end
