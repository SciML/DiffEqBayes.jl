using Distributions
function stan_string(p::Bernoulli,i)
	parameters = (params(p)[1])
	return string("theta[$i] ~ bernoulli($parameters)",";")
end
function stan_string(p::Binomial,i)
	parameters = (params(p)[1])
	return string("theta[$i] ~ binomial($parameters)",";")
end
function stan_string(p::BetaBinomial,i)
	parameters = (params(p)[1],params(p)[2],params(p)[3])
	return string("theta[$i] ~ beta_binomial$parameters",";")
end
function stan_string(p::Hypergeometric,i)
	parameters = (params(p)[1])
	return string("theta[$i] ~ hypergeometric($parameters)",";")
end
function stan_string(p::Categorical,i)
	parameters = (params(p)[1])
	return string("theta[$i] ~ categorical($parameters)",";")
end
function stan_string(p::NegativeBinomial,i)
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ neg_binomial$parameters",";")
end
function stan_string(p::Poisson,i)
	parameters = (params(p)[1])
	return string("theta[$i] ~ beta_binomial($parameters)",";")
end
function stan_string(p::Normal,i)
	parameters = (params(p)[1],params(p)[2])
	string("theta[$i] ~ normal$parameters",";")
end
function stan_string(p::TDist,i) 
	parameters = (params(p)[1])
	return string("theta[$i] ~ normal($parameters,0,1)",";")
end
function stan_string(p::Cauchy,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ cauchy$parameters",";")
end
function stan_string(p::Laplace,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ double_exponential$parameters",";")
end
function stan_string(p::Logistic,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ logistic$parameters",";")
end
function stan_string(p::Gumbel,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ gumbel$parameters",";")
end
function stan_string(p::LogNormal,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ lognormal$parameters",";")
end
function stan_string(p::Chisq,i) 
	parameters = (params(p)[1])
	return string("theta[$i] ~ chi_square($parameters)",";")
end
function stan_string(p::Exponential,i) 
	parameters = (params(p)[1])
	return string("theta[$i] ~ exponential($parameters)",";")
end
function stan_string(p::Gamma,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ gamma$parameters",";")
end
function stan_string(p::InverseGamma,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ inv_gamma$parameters",";")
end
function stan_string(p::Weibull,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ weibull$parameters",";")
end
function stan_string(p::Frechet,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ frechet$parameters",";")
end
function stan_string(p::Rayleigh,i) 
	parameters = (params(p)[1])
	return string("theta[$i] ~ rayleigh($parameters)",";")
end
function stan_string(p::Pareto,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ pareto$parameters",";")
end
function stan_string(p::GeneralizedPareto,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ pareto_type_2$parameters",";")
end
function stan_string(p::Beta,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ beta$parameters",";")
end
function stan_string(p::Uniform,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ uniform$parameters",";")
end
function stan_string(p::VonMises,i) 
	parameters = (params(p)[1],params(p)[2])
	return string("theta[$i] ~ von_mises$parameters",";")
end
function stan_string{T<:Truncated}(p::T,i)
	lower = p.lower
	upper = p.upper
	raw_string = stan_string(p.untruncated,i)
	return string(raw_string[1:end-1]," T[$lower,$upper]", ";")
end