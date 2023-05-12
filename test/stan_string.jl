using Distributions, DiffEqBayes

println("Starting the test")
@test DiffEqBayes.stan_string(Bernoulli(1)) == "bernoulli(1.0)"
@test DiffEqBayes.stan_string(Bernoulli) == "bernoulli"

@test DiffEqBayes.stan_string(Binomial(5, 0.3)) == "binomial(5, 0.3)"
@test DiffEqBayes.stan_string(Binomial) == "binomial"

@test DiffEqBayes.stan_string(BetaBinomial(5, 1, 1)) == "beta_binomial(5, 1.0, 1.0)"
@test DiffEqBayes.stan_string(BetaBinomial) == "beta_binomial"

@test DiffEqBayes.stan_string(Hypergeometric(5, 5, 3)) == "hypergeometric(5, 5, 3)"
@test DiffEqBayes.stan_string(Hypergeometric) == "hypergeometric"

@test DiffEqBayes.stan_string(NegativeBinomial(5, 0.3)) == "neg_binomial(5.0, 0.3)"
@test DiffEqBayes.stan_string(NegativeBinomial) == "neg_binomial"

@test DiffEqBayes.stan_string(Poisson(5)) == "poisson(5.0)"
@test DiffEqBayes.stan_string(Poisson) == "poisson"

@test DiffEqBayes.stan_string(Normal(0, 1)) == "normal(0.0, 1.0)"
@test DiffEqBayes.stan_string(Normal) == "normal"

@test DiffEqBayes.stan_string(TDist(5)) == "student_t(5.0,0,1)"
@test DiffEqBayes.stan_string(TDist) == "student_t"

@test DiffEqBayes.stan_string(Cauchy(0, 1)) == "cauchy(0.0, 1.0)"
@test DiffEqBayes.stan_string(Cauchy) == "cauchy"

@test DiffEqBayes.stan_string(Laplace(0, 1)) == "double_exponential(0.0, 1.0)"
@test DiffEqBayes.stan_string(Laplace) == "double_exponential"

@test DiffEqBayes.stan_string(Distributions.Logistic(0, 1)) == "logistic(0.0, 1.0)"
@test DiffEqBayes.stan_string(Distributions.Logistic) == "logistic"

@test DiffEqBayes.stan_string(Gumbel(0, 1)) == "gumbel(0.0, 1.0)"
@test DiffEqBayes.stan_string(Gumbel) == "gumbel"

@test DiffEqBayes.stan_string(LogNormal(0, 1)) == "lognormal(0.0, 1.0)"
@test DiffEqBayes.stan_string(LogNormal) == "lognormal"

@test DiffEqBayes.stan_string(Chisq(5)) == "chi_square(5.0)"
@test DiffEqBayes.stan_string(Chisq) == "chi_square"

@test DiffEqBayes.stan_string(Exponential(5)) == "exponential(5.0)"
@test DiffEqBayes.stan_string(Exponential) == "exponential"

@test DiffEqBayes.stan_string(Gamma(2, 3)) == "gamma(2.0, 3.0)"
@test DiffEqBayes.stan_string(Gamma) == "gamma"

@test DiffEqBayes.stan_string(InverseGamma(2, 3)) == "inv_gamma(2.0, 3.0)"
@test DiffEqBayes.stan_string(InverseGamma) == "inv_gamma"

@test DiffEqBayes.stan_string(Weibull(1, 1)) == "weibull(1.0, 1.0)"
@test DiffEqBayes.stan_string(Weibull) == "weibull"

@test DiffEqBayes.stan_string(Frechet(1, 1)) == "frechet(1.0, 1.0)"
@test DiffEqBayes.stan_string(Frechet) == "frechet"

@test DiffEqBayes.stan_string(Rayleigh(5)) == "rayleigh(5.0)"
@test DiffEqBayes.stan_string(Rayleigh) == "rayleigh"

@test DiffEqBayes.stan_string(Pareto(2, 3)) == "pareto(2.0, 3.0)"
@test DiffEqBayes.stan_string(Pareto) == "pareto"

@test DiffEqBayes.stan_string(GeneralizedPareto(0, 1, 2)) == "pareto_type_2(0.0, 1.0, 2.0)"
@test DiffEqBayes.stan_string(GeneralizedPareto) == "pareto_type_2"

@test DiffEqBayes.stan_string(Beta(3, 3)) == "beta(3.0, 3.0)"
@test DiffEqBayes.stan_string(Beta) == "beta"

@test DiffEqBayes.stan_string(Uniform(1, 4)) == "uniform(1.0, 4.0)"
@test DiffEqBayes.stan_string(Uniform) == "uniform"

@test DiffEqBayes.stan_string(VonMises(0, 2)) == "von_mises(0.0, 2.0)"
@test DiffEqBayes.stan_string(VonMises) == "von_mises"

@test DiffEqBayes.stan_string(truncated(Normal(1, 2), -1, 3)) ==
      "normal(1.0, 2.0) T[-1.0,3.0]"
@test DiffEqBayes.stan_string(truncated(Normal(1, 2); lower = -1, upper = 3)) ==
      "normal(1.0, 2.0) T[-1.0,3.0]"
@test DiffEqBayes.stan_string(truncated(Normal(1, 2), -Inf, 3)) ==
      "normal(1.0, 2.0) T[,3.0]"
@test DiffEqBayes.stan_string(truncated(Normal(1, 2), nothing, 3)) ==
      "normal(1.0, 2.0) T[,3.0]"
@test DiffEqBayes.stan_string(truncated(Normal(1, 2); upper = 3)) ==
      "normal(1.0, 2.0) T[,3.0]"
@test DiffEqBayes.stan_string(truncated(Normal(1, 2), -1, Inf)) ==
      "normal(1.0, 2.0) T[-1.0,]"
@test DiffEqBayes.stan_string(truncated(Normal(1, 2), -1, nothing)) ==
      "normal(1.0, 2.0) T[-1.0,]"
@test DiffEqBayes.stan_string(truncated(Normal(1, 2); lower = -1)) ==
      "normal(1.0, 2.0) T[-1.0,]"

@test DiffEqBayes.stan_string(truncated(Beta(2, 3), 0.1, 0.4)) ==
      "beta(2.0, 3.0) T[0.1,0.4]"
@test DiffEqBayes.stan_string(truncated(Beta(2, 3); lower = 0.1, upper = 0.4)) ==
      "beta(2.0, 3.0) T[0.1,0.4]"
@test DiffEqBayes.stan_string(truncated(Beta(2, 3), 0, 0.4)) == "beta(2.0, 3.0) T[,0.4]"
@test DiffEqBayes.stan_string(truncated(Beta(2, 3), nothing, 0.4)) ==
      "beta(2.0, 3.0) T[,0.4]"
@test DiffEqBayes.stan_string(truncated(Beta(2, 3); upper = 0.4)) ==
      "beta(2.0, 3.0) T[,0.4]"
@test DiffEqBayes.stan_string(truncated(Beta(2, 3), 0.1, 1)) == "beta(2.0, 3.0) T[0.1,]"
@test DiffEqBayes.stan_string(truncated(Beta(2, 3), 0.1, nothing)) ==
      "beta(2.0, 3.0) T[0.1,]"
@test DiffEqBayes.stan_string(truncated(Beta(2, 3); lower = 0.1)) ==
      "beta(2.0, 3.0) T[0.1,]"
