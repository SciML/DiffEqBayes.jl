function turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors = nothing;
  num_samples=1000, delta=0.65, kwargs...)

  N = length(priors)
  _theta = Vector{Real}(undef, N)
  @model bif(x) = begin
    for i in 1:length(priors)
      _theta[i] ~ priors[i]
    end
    σ ~ InverseGamma(2, 3)

    theta = convert(Array{typeof(first(_theta))}, _theta)
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
      x[:,i] ~ MvNormal(res, σ*ones(length(prob.u0)))
    end
  end

  model = bif(data)
  chn = sample(model, Turing.IS(num_samples))
end