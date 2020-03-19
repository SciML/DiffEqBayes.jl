using Turing: Tracker

function turing_inference(
    prob::DiffEqBase.DEProblem,
    alg,
    t,
    data,
    priors;
    likelihood_dist_priors = [InverseGamma(2, 3)],
    likelihood = (u,p,t,σ) -> MvNormal(u, σ[1]*ones(length(u))),
    num_samples=1000, sampler = Turing.NUTS(0.65),
    syms = [Turing.@varname(theta[i]) for i in 1:length(priors)],
    obsvbls = 1:size(data, 1),
    sample_u0 = false, 
    kwargs...,
)
    N = length(priors)
    Turing.@model mf(x, ::Type{T} = Float64) where {T <: Real} = begin
        theta = Vector{T}(undef, length(priors))
        for i in 1:length(priors)
            theta[i] ~ NamedDist(priors[i], syms[i])
        end
        σ = Vector{T}(undef, length(likelihood_dist_priors))
        for i in 1:length(likelihood_dist_priors)
            σ[i] ~ likelihood_dist_priors[i]
        end
        nu = length(prob.u0)
        u0 = convert.(T, sample_u0 ? theta[1:nu] : prob.u0)
        p = convert.(T, sample_u0 ? theta[(nu + 1):end] : theta)
        _saveat = t === nothing ? Float64[] : t
        p_tmp = remake(prob, u0 = Tracker.data.(u0), p = Tracker.data.(p))
        sol_tmp = solve(p_tmp, alg; saveat = _saveat, kwargs...)
        if T <: Tracker.TrackedReal || Turing.ADBACKEND[] == :zygote
            sol_tmp′ = concrete_solve(prob, alg, u0, p; saveat = _saveat, kwargs...)
        else
            sol_tmp′ = sol_tmp
        end

        if sol_tmp isa DiffEqBase.AbstractEnsembleSolution
            failure = any((s.retcode != :Success for s in sol_tmp)) && any((s.retcode != :Terminated for s in sol_tmp))
        else
            failure = sol_tmp.retcode != :Success && sol_tmp != :Terminated
        end

        if failure
            @logpdf() = -Inf
            return
        end
        if sol_tmp isa DiffEqBase.AbstractNoTimeSolution
            res = sol_tmp′[obsvbls]
            x ~ likelihood(res, theta, Inf, σ)
        else
            for i = 1:length(t)
                res = sol_tmp′[obsvbls, i]
                _t = sol_tmp.t[i]
                x[:, i] ~ likelihood(res, theta, _t, σ)
            end
        end
        return 
    end

    # Instantiate a Model object.
    model = mf(data)
    chn = sample(model, sampler, num_samples)
    return chn
end
