struct DynamicHMCPosterior
    alg
    problem
    likelihood
    priors
    kwargs
end

function (P::DynamicHMCPosterior)(a)
    @unpack alg, problem, likelihood, priors, kwargs = P
    u0_new = Flux.Tracker.TrackedReal.(problem.u0)
    prob = remake(problem,u0 = u0_new,p=a.a)
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
    PT = TransformedLogDensity(transformations, P)
    PTG = FluxGradientLogDensity(PT);

    chain, NUTS_tuned = NUTS_init_tune_mcmc(PTG,num_samples, ϵ=ϵ)
    posterior = transform.(Ref(PTG.transformation), get_position.(chain));

    return posterior, chain, NUTS_tuned
end
