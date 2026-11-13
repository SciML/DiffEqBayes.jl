using DiffEqBase: SDEProblem
using StochasticDiffEq: EM
using Turing: NUTS, arraydist, @model
using Distributions: MvNormal, Poisson, Normal
using LinearAlgebra: Diagonal

"""
    sde_inference(prob::SDEProblem, data, t, priors; kwargs...)

Bayesian parameter estimation for discretely observed SDEs using the
Euler-Maruyama pseudo-likelihood.
"""
function sde_inference(
        prob::SDEProblem, data, t, priors;
        likelihood = :gaussian,
        sigma = 0.1,
        num_samples = 1000,
        dt = 0.1,
        kwargs...
    )
    @model function sde_model(data, t, prob, priors)
        theta ~ arraydist(priors)
        p = theta
        # Euler-Maruyama transition likelihood
        for i in 1:(length(t) - 1)
            delta_t = t[i + 1] - t[i]
            u = data[:, i]
            f_val = similar(u)
            g_val = similar(u)
            prob.f(f_val, u, p, t[i])
            prob.g(g_val, u, p, t[i])
            mu_val = u + f_val * delta_t
            Sigma = Diagonal((g_val .^ 2) .* delta_t .+ sigma^2)
            if likelihood == :gaussian
                data[:, i + 1] ~ MvNormal(mu_val, Sigma)
            elseif likelihood == :poisson
                data[:, i + 1] ~ arraydist(Poisson.(mu_val))
            else
                error("Unsupported likelihood: $likelihood")
            end
        end
    end
    chain = sample(sde_model(data, t, prob, priors), NUTS(), num_samples; kwargs...)
    return chain
end
