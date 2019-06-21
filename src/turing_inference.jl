function turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors = nothing;
  likelihood_dist_priors = [InverseGamma(2, 3)],
  likelihood = (u,p,t,σ) -> MvNormal(u, σ[1]*ones(length(u))),
  num_samples=1000, sampler = Turing.NUTS(num_samples, 0.65),
  kwargs...)

  N = length(priors)
  @model bif(x) = begin
    val ~ priors[1]
    theta = Vector{typeof(val)}(undef,length(priors))
    for i in 1:length(priors)
      theta[i] ~ priors[i]
    end
    val2 ~ likelihood_dist_priors[1]
    σ = Vector{typeof(val)}(undef,length(priors))
    for i in 1:length(likelihood_dist_priors)
      σ[i] ~ likelihood_dist_priors[i]
    end
    p_tmp = remake(prob, u0 = convert.(eltype(theta), (prob.u0)), p = theta)
    sol_tmp = solve(p_tmp, alg; saveat = t, kwargs...)

    fill_length = length(t) - length(sol_tmp.u)
    for i in 1:fill_length
      if eltype(sol_tmp.u) <: Number
        push!(sol_tmp.u, Inf)
      else
        push!(sol_tmp.u, fill(Inf, size(sol_tmp[1])))
      end
    end

    for i = 1:length(t)
      res = sol_tmp.u[i]
      _t   = sol_tmp.t[i]
      x[:,i] ~ likelihood(res,theta,_t,σ)
    end
  end

  model = bif(data)
  chn = sample(model, sampler)
end
