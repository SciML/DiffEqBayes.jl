"""
    sde_inference(prob::SDEProblem, data, t, priors; kwargs...)

Bayesian parameter estimation for discretely observed SDEs using the
Euler-Maruyama pseudo-likelihood.
"""
function sde_inference(
        prob::SDEProblem, data, t, priors;
        likelihood = :gaussian,
        sigma = nothing,
        kwargs...
    )
    @model function sde_model(data, t, prob, priors)
        theta ~ arraydist(priors)
        sol = solve(remake(prob, p = theta), EM(); saveat = t)
        for i in 1:length(t)
            if likelihood == :gaussian
                data[:, i] ~ MvNormal(sol[:, i], sigma)
            elseif likelihood == :poisson
                data[:, i] ~ Poisson.(sol[:, i])
            else
                error("Unsupported likelihood: $likelihood")
            end
        end
    end
    chain = sample(sde_model(data, t, prob, priors), NUTS(), kwargs...)
    return chain
end

