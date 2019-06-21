function turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors = nothing;
  likelihood_dist_priors = [InverseGamma(2, 3)],
  likelihood = (u,p,t,σ) -> MvNormal(u, σ[1]*ones(length(u))),
  num_samples=1000, delta=0.65, kwargs...)

  N = length(priors)
  _theta = Vector{Real}(undef, N)
  _σ = Vector{Real}(undef, length(likelihood_dist_priors))
  @model bif(x) = begin
    for i in 1:length(priors)
      _theta[i] ~ priors[i]
    end
    for i in 1:length(likelihood_dist_priors)
      _σ[i] ~ InverseGamma(2, 3)
    end

    theta = convert(Array{typeof(first(_theta))}, _theta)
    σ = convert(Array{typeof(first(_σ))}, _σ)

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
  chn = sample(model, Turing.IS(num_samples))
end
