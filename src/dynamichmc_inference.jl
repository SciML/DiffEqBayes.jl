struct LotkaVolterraPosterior{Problem, Data, A_Prior, ObservationTimes, ErrorDist}
    problem::Problem
    data::Data
    a_prior::A_Prior
    t::ObservationTimes
    ϵ_dist::ErrorDist
end

function (P::LotkaVolterraPosterior)(θ)
    @unpack problem, data, a_prior, t, ϵ_dist = P
    a = θ[1]
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
    ℓ + logpdf(a_prior, a)
end

function dynamichmc_inference(prob1::DEproblem,data,priors,t,σ = 0.01)
    P = LotkaVolterraPosterior(prob1, data, priors, t, Normal(0.0, σ))

    parameter_transformation = TransformationTuple((bridge(ℝ, ℝ⁺, ))) # assuming a > 0

    PT = TransformLogLikelihood(P, parameter_transformation)
    PTG = ForwardGradientWrapper(PT, zeros(1));

# NOTE: starting from correct parameter is important, otherwise stepsize
# adaptation is not handled well. would probably maximize PT in a real-life
# setting.
    PO = OnceDifferentiable(x -> -P(x), [2.0])
    a₀ = Optim.minimizer(optimize(a -> -P([a]), 0, 10))
    sample, _ = NUTS_init_tune_mcmc(PTG,
                                inverse(parameter_transformation, (a₀, )),
                                1000)

    posterior = ungrouping_map(Vector, get_transformation(PT) ∘ get_position, sample)

    a, = posterior
    a
end