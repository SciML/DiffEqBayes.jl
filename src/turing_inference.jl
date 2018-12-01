function turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors = nothing;
                          num_samples=1000, delta=0.65, kwargs...)

  function bif(vi, sampler, x=data)
    _lp = 0.0
    N = length(priors)
    _theta = Vector(undef,N)

    for i in 1:length(priors)
      _theta[i], __lp = Turing.assume(sampler,
                          priors[i],
                          Turing.VarName(vi, [:bif, Symbol("theta$i")], ""),
                          vi)
      _lp += __lp
    end

    theta = convert(Array{typeof(first(_theta))},_theta)

    σ, __lp = Turing.assume(sampler,
                    InverseGamma(2, 3),
                    Turing.VarName(vi, [:bif, :σ], ""),
                    vi)
    _lp += __lp

    p_tmp = remake(prob, u0=convert.(eltype(theta),(prob.u0)),p=theta); sol_tmp = solve(p_tmp,alg;saveat=t,kwargs...)

    for i = 1:length(t)
      res = sol_tmp.u[i]
      # x[:,i] ~ MvNormal(res, σ*ones(2))
      __lp = Turing.observe(
        sampler,
        MvNormal(res, σ*ones(length(prob.u0))),   # Distribution
        x[:,i],    # Data point
        vi
      )
      _lp += __lp
    end

    vi.logp = _lp
    vi
  end

  bif() = bif(Turing.VarInfo(), nothing)

  chn = sample(bif, Turing.NUTS(num_samples, delta))
end
