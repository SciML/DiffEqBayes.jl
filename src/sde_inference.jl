# Euler-Maruyama transition pseudo-likelihood for discretely observed SDEs.
#
# Between two consecutive observations the SDE `du = f(u, p, t) dt + G(u, p, t) dW` is
# approximated by one Euler-Maruyama step started at the earlier observation, so that
#
#     data[:, i + 1] | data[:, i] ~ N(data[:, i] + f Δt, G Gᵀ Δt + Diagonal(sigma .^ 2))
#
# with `f` and `G` evaluated at `(data[:, i], p, t[i])` and `Δt = t[i + 1] - t[i]`. The
# density is evaluated in closed form, so no trajectory is simulated and the result is a
# deterministic, differentiable function of the parameters.

_em_component(x::Number, _) = x
_em_component(x::AbstractArray, j) = x[j]

# Diagonal noise: `diffusion` holds one coefficient per state (or one shared number).
function _em_transition_logpdf(y, u, drift, diffusion, dt, sigma)
    logp = zero(dt)
    for j in eachindex(y)
        variance = abs2(_em_component(diffusion, j)) * dt + abs2(_em_component(sigma, j))
        variance > 0 || return -Inf
        residual = y[j] - _em_component(u, j) - _em_component(drift, j) * dt
        logp -= (log(2π * variance) + abs2(residual) / variance) / 2
    end
    return logp
end

# Non-diagonal noise: `diffusion` is the `length(u) × m` matrix multiplying `m` Wiener
# processes, so the transition covariance has off-diagonal entries.
function _em_transition_logpdf(y, u, drift, diffusion::AbstractMatrix, dt, sigma)
    noise_matrix = Matrix(diffusion)
    nugget = [abs2(_em_component(sigma, j)) for j in axes(noise_matrix, 1)]
    covariance = (noise_matrix * transpose(noise_matrix)) .* dt + Diagonal(nugget)
    all(isfinite, covariance) || return -Inf
    factorization = cholesky(Symmetric(covariance); check = false)
    issuccess(factorization) || return -Inf
    residual = y .- u .- drift .* dt
    mahalanobis = dot(residual, factorization \ residual)
    return -(length(residual) * log(2π) + logdet(factorization) + mahalanobis) / 2
end

# Sum of the transition log-densities over all consecutive observation pairs. `T` is the
# element type of the sampled parameters; it sizes the in-place buffers so that automatic
# differentiation can write dual numbers into them.
function _em_loglikelihood(
        prob::SciMLBase.SDEProblem, p, t, data, sigma, ::Type{T}
    ) where {T}
    S = promote_type(T, float(eltype(data)))
    logp = zero(S)
    if SciMLBase.isinplace(prob)
        prototype = prob.noise_rate_prototype
        drift = Vector{S}(undef, size(data, 1))
        diffusion = if prototype === nothing
            Vector{S}(undef, size(data, 1))
        else
            fill!(similar(prototype, S), zero(S))
        end
        for i in 1:(length(t) - 1)
            u = data[:, i]
            prob.f(drift, u, p, t[i])
            prob.g(diffusion, u, p, t[i])
            logp += _em_transition_logpdf(
                view(data, :, i + 1), u, drift, diffusion, t[i + 1] - t[i], sigma
            )
        end
    else
        scalar_state = prob.u0 isa Number
        for i in 1:(length(t) - 1)
            u = scalar_state ? data[1, i] : data[:, i]
            drift = prob.f(u, p, t[i])
            diffusion = prob.g(u, p, t[i])
            logp += _em_transition_logpdf(
                view(data, :, i + 1), u, drift, diffusion, t[i + 1] - t[i], sigma
            )
        end
    end
    return logp
end

function _em_loglikelihood(prob::SciMLBase.SDEProblem, p, t, data, sigma)
    return _em_loglikelihood(prob, p, t, data, sigma, eltype(p))
end

function _check_sde_inference_inputs(prob::SciMLBase.SDEProblem, t, data, sigma)
    u0 = prob.u0
    if !(u0 isa Union{Number, AbstractVector})
        throw(ArgumentError("`sde_inference` requires a scalar or vector state `u0`."))
    end
    n = length(u0)
    prototype = prob.noise_rate_prototype
    if prob.noise !== nothing
        throw(ArgumentError("`sde_inference` does not support custom `noise` processes."))
    elseif length(t) < 2
        throw(ArgumentError("`sde_inference` needs at least two observation times."))
    elseif !all(i -> t[i + 1] > t[i], 1:(length(t) - 1))
        throw(ArgumentError("The observation times `t` must be strictly increasing."))
    elseif size(data) != (n, length(t))
        expected = (n, length(t))
        throw(DimensionMismatch("`data` must have size $expected, got $(size(data))."))
    elseif !(sigma isa Number || length(sigma) == n)
        throw(DimensionMismatch("`sigma` must be a number or have $n entries."))
    elseif !all(s -> isfinite(s) && s >= 0, sigma)
        throw(ArgumentError("`sigma` must be finite and non-negative."))
    elseif prototype !== nothing && size(prototype, 2) < n && all(iszero, sigma)
        # Fewer Wiener processes than states: `G Gᵀ Δt` is singular for every parameter.
        message = "The problem has fewer noise terms than states, so the " *
            "Euler-Maruyama transition covariance is singular. Pass `sigma > 0`."
        throw(ArgumentError(message))
    end
    return nothing
end

"""
    sde_inference(prob::SDEProblem, t, data, priors; kwargs...)

Run Bayesian parameter inference for a discretely observed stochastic differential
equation with Turing.jl, using the Euler-Maruyama pseudo-likelihood.

For `du = f(u, p, t) dt + G(u, p, t) dW`, each transition between consecutive
observations is approximated by one Euler-Maruyama step started at the earlier
observation:

```math
x_{i+1} \\mid x_i \\sim \\mathcal{N}\\left(x_i + f(x_i, p, t_i) \\Delta t_i,\\;
G(x_i, p, t_i) G(x_i, p, t_i)^\\top \\Delta t_i + \\operatorname{diag}(\\sigma^2)\\right),
```

where `xᵢ = data[:, i]` and `Δtᵢ = t[i + 1] - t[i]`. The density is evaluated in closed
form from `prob.f` and `prob.g`. The SDE is never solved, so no solver, `dt`, or
`StochasticDiffEq` is needed, and the log-density is a deterministic function of the
parameters that gradient-based samplers such as NUTS can target.

The approximation error grows with the spacing of the observations. It is accurate when
`Δtᵢ` is small relative to the time scales of the drift and diffusion, and it biases the
posterior when the data are sparse in time.

# Arguments

- `prob`: an `SDEProblem`. Only its drift, diffusion, and `noise_rate_prototype` are
  used; `u0` only sets the shape of the state (a number or a vector), and `tspan` is
  ignored. In-place and out-of-place problems, diagonal noise, and non-diagonal noise
  (through `noise_rate_prototype`) are supported. Custom `noise` processes are not.
- `t`: strictly increasing observation times.
- `data`: observations, with `data[:, i]` the full state observed at `t[i]`. A vector is
  accepted for a one-dimensional state. The first observation is conditioned on.
- `priors`: an iterable of priors, one per parameter. The sampled values are passed to
  the drift and diffusion as the parameter vector `p`; a structured `prob.p` is rebuilt
  from them through `SciMLStructures.Tunable`.

# Keywords

- `sigma`: standard deviation of an independent Gaussian term added to every transition,
  as a number or one value per state. The default `0` gives the plain Euler-Maruyama
  transition density. A positive value regularizes rank-deficient diffusions (fewer
  noise terms than states, or noise-free components) and crudely accounts for
  measurement error. It is not a state-space model: the previous observation is still
  treated as the exact state.
- `syms`: Turing variable names for the sampled parameters. The default creates one
  name per entry in `priors`.
- `progress`: whether Turing displays sampling progress.
- `sample_args`: named arguments controlling Turing's sampler, with defaults of
  `NUTS(0.65)`, `MCMCSerial()`, `1000` samples, and one chain.
- `sample_kwargs`: additional keyword arguments forwarded to `Turing.sample`.

# Returns

Returns the chain produced by `Turing.sample`, indexed using the `VarName` values in
`syms`, as for [`turing_inference`](@ref).

# Examples

```julia
drift(u, p, t) = p[1] .* (p[2] .- u)
diffusion(u, p, t) = [p[3]]
prob = SDEProblem(drift, diffusion, [0.0], (0.0, 1.0), [1.0, 0.0, 0.5])
priors = [Uniform(0, 5), Normal(0, 1), Uniform(0, 2)]
chain = sde_inference(prob, t, data, priors; syms = [:rate, :level, :noise])
mean(chain[@varname(rate)])
```
"""
function sde_inference(
        prob::SciMLBase.SDEProblem,
        t,
        data,
        priors;
        sigma = 0.0,
        syms = [Turing.@varname(theta[i]) for i in 1:length(priors)],
        progress = false,
        sample_args = NamedTuple(), # accept Turing.jl sample args
        sample_kwargs = Dict() # accept Turing.jl sample kwargs
    )
    observations = data isa AbstractVector ? reshape(data, 1, :) : data
    _check_sde_inference_inputs(prob, t, observations, sigma)
    param_vns = map(_as_varname, syms)
    # default args are updated with user supplied args
    sample_args = (;
        sampler = Turing.NUTS(0.65),
        parallel_type = MCMCSerial(),
        num_samples = 1000,
        n_chains = 1,
        sample_args...,
    )

    # Plain arrays (and unset parameters) receive the sampled vector directly; structured
    # containers are rebuilt from it, as in `turing_inference`.
    repack = if prob.p isa AbstractArray || !SciMLStructures.isscimlstructure(prob.p)
        identity
    else
        SciMLStructures.canonicalize(SciMLStructures.Tunable(), prob.p)[2]
    end

    Turing.@model function sde_infer(x, ::Type{T} = Float64) where {T <: Real}
        theta = Vector{T}(undef, length(priors))
        for i in 1:length(priors)
            # Same run-time varname handling as in `turing_inference`.
            param, __varinfo__ = tilde_assume!!(
                __model__.context, priors[i], param_vns[i], NoTemplate(), __varinfo__
            )
            theta[i] = param
        end
        Turing.@addlogprob! _em_loglikelihood(prob, repack(theta), t, x, sigma, T)
        return
    end false

    chn = sample(
        sde_infer(observations),
        sample_args...;
        progress = progress,
        sample_kwargs...
    )
    return chn
end
