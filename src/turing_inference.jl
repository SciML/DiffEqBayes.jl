function turing_inference(
        prob::DiffEqBase.DEProblem,
        alg,
        t,
        data,
        priors;
        likelihood_dist_priors = [InverseGamma(2, 3)],
        likelihood = (u, p, t, σ) -> MvNormal(
            u, σ[1]^2 * Diagonal(ones(length(u)))
        ),
        syms = [Turing.@varname(theta[i]) for i in 1:length(priors)],
        sample_u0 = false,
        progress = false,
        solve_kwargs = Dict(), # accept SciML DiffEq solve kwargs
        sample_args = NamedTuple(), # accept Turing.jl sample args
        sample_kwargs = Dict() # accept Turing.jl sample kwargs
)
    N = length(priors)
    # default args are updated with user supplied args
    solve_kwargs = merge(Dict(:save_idxs => nothing), solve_kwargs)
    sample_args = (;
        sampler = Turing.NUTS(0.65),
        parallel_type = MCMCSerial(),
        num_samples = 1000,
        n_chains = 1,
        sample_args...
    )

    _p, repack, aliases = if SciMLStructures.isscimlstructure(prob.p)
        SciMLStructures.canonicalize(SciMLStructures.Tunable(), prob.p)
    else
        prob.p, identity, true
    end

    Turing.@model function infer(x, ::Type{T} = Float64) where {T <: Real}
        theta = Vector{T}(undef, length(priors))
        for i in 1:length(priors)
            theta[i] ~ NamedDist(priors[i], syms[i])
        end
        σ = Vector{T}(undef, length(likelihood_dist_priors))
        for i in 1:length(likelihood_dist_priors)
            σ[i] ~ likelihood_dist_priors[i]
        end
        nu = solve_kwargs[:save_idxs] === nothing ? length(prob.u0) :
             length(solve_kwargs[:save_idxs])
        u0 = convert.(T, sample_u0 ? theta[1:nu] : prob.u0)
        p = convert.(T, sample_u0 ? theta[(nu + 1):end] : theta)
        if length(u0) < length(prob.u0)
            # assumes u is ordered such that the observed variables are in the beginning, consistent with ordered theta
            for i in length(u0):length(prob.u0)
                push!(u0, convert(T, prob.u0[i]))
            end
        end
        _saveat = t === nothing ? Float64[] : t
        sol = solve(prob, alg; u0 = u0, p = repack(p), saveat = _saveat,
            progress = progress, solve_kwargs...)
        failure = size(sol, 2) < length(_saveat)

        if failure
            Turing.DynamicPPL.acclogp!!(__varinfo__, -Inf)
            return
        end
        if ndims(sol) == 1
            x ~ likelihood(Array(sol), theta, Inf, σ)
        else
            for i in 1:length(t)
                x[:, i] ~ likelihood(sol[:, i], theta, sol.t[i], σ)
            end
        end
        return
    end false

    # Instantiate a Model object.
    model = infer(data)
    chn = sample(
        model,
        sample_args...;
        progress = progress,
        sample_kwargs...
    )
    return chn
end
