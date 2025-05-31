"""
$(TYPEDEF)
Defines a callable that returns the log density for given parameter values when called with
a `NamedTuple` `(parameters = ..., σ = ...)` where `parameters` is a vector of parameters,
and `σ` is the vector of noise scales.
For a common use case, see [`dynamichmc_inference`](@ref).
# Fields
$(FIELDS)
"""
Base.@kwdef struct DynamicHMCPosterior{TA, TP, TD, TT, TR, TS, TK, TI, TRe}
    "Algorithm for the ODE solver."
    algorithm::TA
    "An ODE problem definition (`DiffEqBase.DEProblem`)."
    problem::TP
    "Time values at which the simulated path is compared to `data`."
    t::TT
    "Data, as a matrix with each time value in a column."
    data::TD
    "Priors for parameters, an iterable with the same length as the number of parameters."
    parameter_priors::TR
    """
    Priors for the noise scales (currently the standard deviation of a normal distribution),
    one for each variable.
    """
    σ_priors::TS
    "Keyword arguments passed on the the ODE solver `solve`."
    solve_kwargs::TK
    sample_u0::Bool
    save_idxs::TI
    repack::TRe
end

function (P::DynamicHMCPosterior)(θ)
    @unpack parameters, σ = θ
    @unpack algorithm, problem, data, t, parameter_priors = P
    @unpack σ_priors, solve_kwargs, sample_u0, save_idxs = P
    T = eltype(parameters)
    nu = save_idxs === nothing ? length(problem.u0) : length(save_idxs)
    u0 = convert.(T, sample_u0 ? parameters[1:nu] : problem.u0)
    p = convert.(T, sample_u0 ? parameters[(nu + 1):end] : parameters)
    if length(u0) < length(problem.u0)
        # assumes u is ordered such that the observed variables are in the beginning, consistent with ordered theta
        for i in length(u0):length(problem.u0)
            push!(u0, convert(T, problem.u0[i]))
        end
    end
    _saveat = t === nothing ? Float64[] : t
    sol = solve(
        problem, algorithm; u0 = u0, p = P.repack(p), saveat = _saveat, save_idxs = save_idxs,
        solve_kwargs...)
    failure = size(sol, 2) < length(_saveat)
    failure && return T(0) * sum(σ) + T(-Inf)
    log_likelihood = sum(sum(map(logpdf, Normal.(0.0, σ), sol[:, i] .- data[:, i]))
    for (i, t) in enumerate(t))
    log_prior_parameters = sum(map(logpdf, parameter_priors, parameters))
    log_prior_σ = sum(map(logpdf, σ_priors, σ))
    log_likelihood + log_prior_parameters + log_prior_σ
end

# function (P::DynamicHMCPosterior)(θ)
#   @unpack parameters, σ = θ
#   @unpack algorithm, problem, data, t, parameter_priors, σ_priors, solve_kwargs = P
#   prob = remake(problem, u0 = convert.(eltype(parameters), problem.u0), p = parameters)
#   solution = solve(prob, algorithm; solve_kwargs...)
#   any((s.retcode ≠ :Success && s.retcode ≠ :Terminated) for s in solution) && return -Inf
#   log_likelihood = sum(sum(logpdf.(Normal.(0.0, σ), solution(t) .- data[:, i]))
#                        for (i, t) in enumerate(t))
#   log_prior_parameters = sum(map(logpdf, parameter_priors, parameters))
#   log_prior_σ = sum(map(logpdf, σ_priors, σ))
#   log_likelihood + log_prior_parameters + log_prior_σ
# end

"""
$(SIGNATURES)
Run MCMC for an ODE problem. Return a `NamedTuple`, which is similar to the one returned by
`DynamicHMC.mcmc_with_warmup`, with an added field `posterior` which contains a vector of
posterior values (transformed from `ℝⁿ`).
# Arguments
- `problem` is the ODE problem
- `algorithm` is the ODE algorithm
- `t` is the time values at which the solution is compared to `data`
- `data` is a matrix of data, with one column for each element in `t`
- `parameter_priors` is an iterable with the length of the number of paramers, and is used
  as a prior on it, should have comparable structure.
- `parameter_transformations`: a `TransformVariables` transformation to mapping `ℝⁿ` to the
  vector of valid parameters.
# Keyword arguments
- `rng` is the random number generator used for MCMC. Defaults to the global one.
- `num_samples` is the number of MCMC draws (default: 1000)
- `AD_gradient_kind` is passed on to `LogDensityProblems.ADgradient`, make sure to `import`
  the corresponding library.
- `solve_kwargs` is passed on to `solve`
- `mcmc_kwargs` are passed on as keyword arguments to `DynamicHMC.mcmc_with_warmup`
"""
function dynamichmc_inference(problem::DiffEqBase.DEProblem, algorithm, t, data,
        parameter_priors,
        parameter_transformations = as(Vector, asℝ₊,
            length(parameter_priors));
        σ_priors = fill(Normal(0, 5), size(data, 1)),
        sample_u0 = false, rng = Random.GLOBAL_RNG,
        num_samples = 1000, AD_gradient_kind = Val(:ForwardDiff),
        save_idxs = nothing, solve_kwargs = (),
        mcmc_kwargs = (initialization = (q = zeros(length(parameter_priors) +
                                           (save_idxs ===
                                            nothing ?
                                            length(data[:, 1]) :
                                            length(save_idxs))),),))
    
        _p, repack, aliases = if SciMLStructures.isscimlstructure(problem.p)
            SciMLStructures.canonicalize(SciMLStructures.Tunable(), problem.p)
        else
            problem.p, identity, true
        end
    
        P = DynamicHMCPosterior(; algorithm = algorithm, problem = problem, t = t, data = data,
        parameter_priors = parameter_priors, σ_priors = σ_priors,
        solve_kwargs = solve_kwargs, sample_u0 = sample_u0,
        save_idxs = save_idxs, repack = repack)


    trans = as((parameters = parameter_transformations,
        σ = as(Vector, asℝ₊, length(σ_priors))))
    ℓ = TransformedLogDensity(trans, P)
    ∇ℓ = LogDensityProblemsAD.ADgradient(AD_gradient_kind, ℓ)
    results = mcmc_with_warmup(rng, ∇ℓ, num_samples; mcmc_kwargs...)
    chain = if haskey(results, :chain) # DynamicHMC < 3.3.0
        results.chain
    else
        eachcol(results.posterior_matrix)
    end
    posterior = map(Base.Fix1(TransformVariables.transform, trans), chain)
    merge((; posterior), results)
end
