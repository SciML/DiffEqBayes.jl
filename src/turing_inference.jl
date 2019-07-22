function turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors;
  likelihood_dist_priors = [InverseGamma(2, 3)],
  likelihood = (u,p,t,σ) -> MvNormal(u, σ[1]*ones(length(u))),
  num_samples=1000, sampler = Turing.NUTS(num_samples, 0.65),
  syms = [Symbol("theta$i") for i in 1:length(priors)],
  kwargs...)

  N = length(priors)
  mf(vi,sampler,model) = begin
    # If x is provided, use the provided values.
    # Otherwise, treat x as an empty vector with
    # two entries.
    x = model.defaults.x

    N = length(priors)
    _theta = Vector(undef,length(priors))

    for i in 1:length(priors)
      _theta[i], __lp = Turing.assume(sampler,
                        priors[i],
                        Turing.VarName([:mf, syms[i]], ""),
                        vi)
      vi.logp += __lp
    end

    theta = convert(Array{typeof(first(_theta))},_theta)
    _σ = Vector(undef,length(likelihood_dist_priors))

    for i in 1:length(likelihood_dist_priors)
      _σ[i], __lp = Turing.assume(sampler,
                      likelihood_dist_priors[i],
                      Turing.VarName([:mf, Symbol("σ$i")], ""),
                      vi)
      vi.logp += __lp
    end
    σ = convert(Array{typeof(first(_σ))},_σ)

    p_tmp = remake(prob, u0 = convert.(eltype(theta), (prob.u0)), p = theta)
    _saveat = t === nothing ? Float64[] : t
    sol_tmp = solve(p_tmp, alg; saveat = _saveat, kwargs...)

    if any((s.retcode != :Success for s in sol_tmp)) && any((s.retcode != :Terminated for s in sol_tmp))
      vi.logp -= Inf
    else
      if sol_tmp isa DiffEqBase.AbstractNoTimeSolution
        res = sol_tmp.u
        __lp = Turing.observe(
          sampler,
          likelihood(res,theta,Inf,σ),   # Distribution
          x,    # Data point
          vi
        )
        vi.logp += __lp
      else
        for i = 1:length(t)
          res = sol_tmp.u[i]
          _t   = sol_tmp.t[i]
          __lp = Turing.observe(
            sampler,
            likelihood(res,theta,_t,σ),   # Distribution
            x[:,i],    # Data point
            vi
          )
          vi.logp += __lp
        end
      end
    end
    vi
  end

  bigtup = Tuple{([[syms[i] for i in 1:length(priors)]; [Symbol("σ$i") for i in 1:length(likelihood_dist_priors)]])...}

  # Define the default value for x when missing
  defaults = (x = data,)

  # Instantiate a Model object.
  model = Turing.Model{bigtup, Tuple{:x}}(mf, data, defaults)

  chn = sample(model, sampler)
end
