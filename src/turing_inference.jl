function turing_inference(prob::DEProblem,alg,t,data,priors = nothing;
                            num_samples=1000, epsilon = 0.02, tau = 4, kwargs...)

  bif(vi, sampler, x=data) = begin
    N = length(priors)
    _theta = Vector(N)
    for i in 1:length(priors)
      _theta[i] = Turing.assume(sampler,
                          priors[i],
                          Turing.VarName(vi, [:bif, Symbol("theta$i")], ""),
                          vi)
    end
    theta = convert(Array{typeof(first(_theta))},_theta)
    σ = Turing.assume(sampler,
                    InverseGamma(2, 3),
                    Turing.VarName(vi, [:bif, :σ], ""),
                    vi)
    p_tmp = problem_new_parameters(prob, theta); sol_tmp = solve(p_tmp,alg;kwargs...)

    for i = 1:length(t)
      res = sol_tmp(t[i])
      # x[:,i] ~ MvNormal(res, σ*ones(2))
      Turing.observe(
        sampler,
        MvNormal(res, σ*ones(2)),   # Distribution
        x[:,i],    # Data point
        vi
      )
    end
    vi
  end

  bif() = bif(Turing.VarInfo(), nothing)

  chn = sample(bif, Turing.HMC(num_samples, epsilon, tau))
end
