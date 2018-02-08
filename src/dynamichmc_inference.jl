struct DynamicHMCPosterior
    problem
    data
    a_prior
    t
    ϵ_dist
end

function (P::DynamicHMCPosterior)(θ)
    #println(θ)
    @unpack problem, data, a_prior, t, ϵ_dist = P
    a = θ
    try
        prob = problem_new_parameters(problem, a)
        sol = solve(prob, Tsit5())
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

function dynamichmc_inference(prob::DEProblem,data,priors,t,transformations,initial,σ = 0.01)
    P = DynamicHMCPosterior(prob, data, priors, t, Normal(0.0, σ))
    
    transformations_tuple = Tuple(transformations)
    parameter_transformation = TransformationTuple(transformations_tuple) # assuming a > 0
    PT = TransformLogLikelihood(P, parameter_transformation)
    PTG = ForwardGradientWrapper(PT, zeros(length(priors)));

    # NOTE: starting from correct parameter is important, otherwise stepsize
    # adaptation is not handled well. would probably maximize PT in a real-life
    # setting
    lower_bound = Float64[]
    upper_bound = Float64[]
    for i in priors
        push!(lower_bound,minimum(i))
        push!(upper_bound,maximum(i))
    end
    optimized = Optim.minimizer(optimize(a -> -P(a),initial,lower_bound,upper_bound,Fminbox{GradientDescent}()))
    #inverse_transforms = Float64[]
    #for i in 1:length(initial)
    #    para = TransformationTuple(transformations[i])
    #    push!(inverse_transforms,inverse(para, (optimized[i], ))[1])
    #end
    #println(inverse_transforms)
    sample, _ = NUTS_init_tune_mcmc(PTG,
                                optimized,
                                1000)
    
    posterior = ungrouping_map(Vector, get_transformation(PT) ∘ get_position, sample)

    posterior
end