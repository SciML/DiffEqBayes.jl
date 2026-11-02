using DiffEqBase: SDEProblem
using StochasticDiffEq: EM
using Turing: NUTS, arraydist

"""
    sde_inference(prob::SDEProblem, data, t, priors; kwargs...)

Bayesian parameter estimation for discretely observed SDEs using the
Euler-Maruyama pseudo-likelihood. This is the SDE branch that complements
the ODE-only `turing_inference` in DiffEqBayes.jl.

# Arguments
- `prob`: An `SDEProblem`.
- `data`: Observed data (array of observations).
- `t`: Observation times.
- `priors`: Prior distributions for the parameters.
- `kwargs`: Additional arguments passed to `Turing.sample`.
"""
function sde_inference(prob::SDEProblem, data, t, priors;
        likelihood = :gaussian,
        sigma = nothing,
        kwargs...)
    # Define the Turing model using the Euler-Maruyama pseudo-likelihood
    @model function sde_model(data, t, prob, priors)
        # Sample parameters from priors
        theta ~ arraydist(priors)

        # Solve the SDE with the current parameters
        sol = solve(remake(prob, p = theta), EM(); saveat = t)

        # Compute the log-likelihood
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

    # Sample using NUTS
    chain = sample(sde_model(data, t, prob, priors), NUTS(), kwargs...)
    return chain
end