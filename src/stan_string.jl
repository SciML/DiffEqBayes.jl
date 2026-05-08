using Distributions

# `stan_string(::Type{D})` returns the bare distribution name (no parameters)
# while `stan_string(::D)` formats the realised parameters. The two are split
# because `Distributions.params` is only defined on instances, so accepting
# `Union{Type{D}, D}` in a single method produces a method-error path that
# JET (correctly) flags as unreachable for the `Type{D}` branch.

function stan_string(::Type{Bernoulli})
    return "bernoulli"
end
function stan_string(p::Bernoulli)
    parameters = (params(p)[1])
    return string("bernoulli($parameters)")
end

function stan_string(::Type{Binomial})
    return "binomial"
end
function stan_string(p::Binomial)
    parameters = (params(p)[1], params(p)[2])
    return string("binomial$parameters")
end

function stan_string(::Type{BetaBinomial})
    return "beta_binomial"
end
function stan_string(p::BetaBinomial)
    parameters = (params(p)[1], params(p)[2], params(p)[3])
    return string("beta_binomial$parameters")
end

function stan_string(::Type{Hypergeometric})
    return "hypergeometric"
end
function stan_string(p::Hypergeometric)
    parameters = (params(p)[1], params(p)[2], params(p)[3])
    return string("hypergeometric$parameters")
end

function stan_string(::Type{<:Categorical})
    return "categorical"
end
function stan_string(p::Categorical)
    parameters = (params(p)[1])
    return string("categorical($parameters)")
end

function stan_string(::Type{NegativeBinomial})
    return "neg_binomial"
end
function stan_string(p::NegativeBinomial)
    parameters = (params(p)[1], params(p)[2])
    return string("neg_binomial$parameters")
end

function stan_string(::Type{Poisson})
    return "poisson"
end
function stan_string(p::Poisson)
    parameters = (params(p)[1])
    return string("poisson($parameters)")
end

function stan_string(::Type{Normal})
    return "normal"
end
function stan_string(p::Normal)
    parameters = (params(p)[1], params(p)[2])
    return string("normal$parameters")
end

function stan_string(::Type{TDist})
    return "student_t"
end
function stan_string(p::TDist)
    parameters = (params(p)[1])
    return string("student_t($parameters,0,1)")
end

function stan_string(::Type{Cauchy})
    return "cauchy"
end
function stan_string(p::Cauchy)
    parameters = (params(p)[1], params(p)[2])
    return string("cauchy$parameters")
end

function stan_string(::Type{Laplace})
    return "double_exponential"
end
function stan_string(p::Laplace)
    parameters = (params(p)[1], params(p)[2])
    return string("double_exponential$parameters")
end

function stan_string(::Type{Distributions.Logistic})
    return "logistic"
end
function stan_string(p::Distributions.Logistic)
    parameters = (params(p)[1], params(p)[2])
    return string("logistic$parameters")
end

function stan_string(::Type{Gumbel})
    return "gumbel"
end
function stan_string(p::Gumbel)
    parameters = (params(p)[1], params(p)[2])
    return string("gumbel$parameters")
end

function stan_string(::Type{LogNormal})
    return "lognormal"
end
function stan_string(p::LogNormal)
    parameters = (params(p)[1], params(p)[2])
    return string("lognormal$parameters")
end

function stan_string(::Type{Chisq})
    return "chi_square"
end
function stan_string(p::Chisq)
    parameters = (params(p)[1])
    return string("chi_square($parameters)")
end

function stan_string(::Type{Exponential})
    return "exponential"
end
function stan_string(p::Exponential)
    parameters = (params(p)[1])
    return string("exponential($parameters)")
end

function stan_string(::Type{Gamma})
    return "gamma"
end
function stan_string(p::Gamma)
    parameters = (params(p)[1], params(p)[2])
    return string("gamma$parameters")
end

function stan_string(::Type{InverseGamma})
    return "inv_gamma"
end
function stan_string(p::InverseGamma)
    parameters = (params(p)[1], params(p)[2])
    return string("inv_gamma$parameters")
end

function stan_string(::Type{Weibull})
    return "weibull"
end
function stan_string(p::Weibull)
    parameters = (params(p)[1], params(p)[2])
    return string("weibull$parameters")
end

function stan_string(::Type{Frechet})
    return "frechet"
end
function stan_string(p::Frechet)
    parameters = (params(p)[1], params(p)[2])
    return string("frechet$parameters")
end

function stan_string(::Type{Rayleigh})
    return "rayleigh"
end
function stan_string(p::Rayleigh)
    parameters = (params(p)[1])
    return string("rayleigh($parameters)")
end

function stan_string(::Type{Pareto})
    return "pareto"
end
function stan_string(p::Pareto)
    parameters = (params(p)[1], params(p)[2])
    return string("pareto$parameters")
end

function stan_string(::Type{GeneralizedPareto})
    return "pareto_type_2"
end
function stan_string(p::GeneralizedPareto)
    parameters = (params(p)[1], params(p)[2], params(p)[3])
    return string("pareto_type_2$parameters")
end

function stan_string(::Type{Beta})
    return "beta"
end
function stan_string(p::Beta)
    parameters = (params(p)[1], params(p)[2])
    return string("beta$parameters")
end

function stan_string(::Type{Uniform})
    return "uniform"
end
function stan_string(p::Uniform)
    parameters = (params(p)[1], params(p)[2])
    return string("uniform$parameters")
end

function stan_string(::Type{VonMises})
    return "von_mises"
end
function stan_string(p::VonMises)
    parameters = (params(p)[1], params(p)[2])
    return string("von_mises$parameters")
end

function stan_string(p::Truncated)
    min_truncated, max_truncated = extrema(p)
    min_untruncated, max_untruncated = extrema(p.untruncated)
    lower = min_truncated == min_untruncated ? "" : string(min_truncated)
    upper = max_truncated == max_untruncated ? "" : string(max_truncated)
    raw_string = stan_string(p.untruncated)
    return string(raw_string, " T[", lower, ",", upper, "]")
end
