struct StanModel{R,C}
  return_code::R
  chain_results::C
end

struct StanODEData
end

function generate_differential_equation(f)
  theta_ex = MacroTools.postwalk(f.fex) do x
    if typeof(x) <: Expr && x.args[1] == :internal_var___p
      return Symbol("theta[$(x.args[2])]")
    else
      return x
    end
  end
  differential_equation = ""
  for i in 1:length(theta_ex.args)-1
    differential_equation = string(differential_equation,theta_ex.args[i], ";\n")
  end
  return differential_equation
end

function generate_priors(f,priors)
  priors_string = ""
  params = f.params
  if priors==nothing
    for i in 1:length(params)
      priors_string = string(priors_string,"theta[$i] ~ normal(0, 1)", " ; ")
    end
  else
    for i in 1:length(params)
      priors_string = string(priors_string,"theta[$i] ~",stan_string(priors[i]),";")
    end
  end
  priors_string
end

function generate_theta(n,priors)
  theta = ""
  for i in 1:n
    upper_bound = ""
    lower_bound = ""
    if maximum(priors[i]) != Inf
      upper_bound = string("upper=",maximum(priors[i]))
    end
    if minimum(priors[i]) != -Inf
      lower_bound = string("lower=",minimum(priors[i]))
    end
    if lower_bound != "" && upper_bound != ""
      theta = string(theta,"real","<$lower_bound",",","$upper_bound>"," theta$i",";")
    elseif lower_bound != ""
      theta = string(theta,"real","<$lower_bound",">"," theta$i",";")
    elseif upper_bound != ""
      theta = string(theta,"real","<","$upper_bound>"," theta$i",";")
    else
      theta = string(theta,"real"," theta$i",";")
    end
  end
  return theta
end

function stan_inference(prob::DEProblem,t,data,priors = nothing;alg=:rk45,
                            num_samples=1000, num_warmup=1000, reltol=1e-3,
                            abstol=1e-6, maxiter=Int(1e5),likelihood=Normal,
                            vars=(StanODEData(),InverseGamma(3,3)))
  length_of_y = length(prob.u0)
  length_of_params = length(vars)
  f = prob.f
  length_of_parameter = string(length(f.params))
  if alg ==:rk45
    algorithm = "integrate_ode_rk45"
  elseif alg == :bdf
    algorithm = "integrate_ode_bdf"
  else
    error("The choices for alg are :rk45 or :bdf")
  end
  hyper_params = ""
  tuple_hyper_params = ""
  setup_params = ""
  thetas = ""
  theta_string = generate_theta(length(priors),priors)
  for i in 1:length(priors)
    thetas = string(thetas,"theta[$i] <- theta$i",";")
  end
  for i in 1:length_of_params
    if isa(vars[i],StanODEData)
      tuple_hyper_params = string(tuple_hyper_params,"u_hat[t,:]",",")
    else
      dist = stan_string(vars[i])
      hyper_params = string(hyper_params,"sigma$(i-1) ~ $dist;")
      tuple_hyper_params = string(tuple_hyper_params,"sigma$(i-1)",",")
      setup_params = string(setup_params,"row_vector<lower=0>[$length_of_y] sigma$(i-1);")
    end
  end
  tuple_hyper_params = tuple_hyper_params[1:endof(tuple_hyper_params)-1]
  differential_equation = generate_differential_equation(f)
  priors_string = string(generate_priors(f,priors))
  stan_likelihood = stan_string(likelihood)
  parameter_estimation_model = "
  functions {
    real[] sho(real t,real[] internal_var___u,real[] theta,real[] x_r,int[] x_i) {
      real internal_var___du[$length_of_y];
      $differential_equation
      return internal_var___du;
      }
    }
  data {
    real u0[$length_of_y];
    int<lower=1> T;
    real internal_var___u[T,$length_of_y];
    real t0;
    real ts[T];
  }
  transformed data {
    real x_r[0];
    int x_i[0];
  }
  parameters {
    $setup_params
    $theta_string
  }
  transformed parameters{
    real theta[$length_of_parameter];
    $thetas
  }
  model{
    real u_hat[T,$length_of_y];
    $hyper_params
    $priors_string
    u_hat = $algorithm(sho, u0, t0, ts, theta, x_r, x_i, $reltol, $abstol, $maxiter);
    for (t in 1:T){
      internal_var___u[t,:] ~ $stan_likelihood($tuple_hyper_params);
      }
  }
  "
  stanmodel = Stanmodel(num_samples=num_samples, num_warmup=num_warmup, name="parameter_estimation_model", model=parameter_estimation_model);
  parameter_estimation_data = Dict("u0"=>prob.u0, "T" => length(t), "internal_var___u" => data', "t0" => prob.tspan[1], "ts" => t)
  return_code, chain_results = stan(stanmodel, [parameter_estimation_data]; CmdStanDir=CMDSTAN_HOME)
  return StanModel(return_code,chain_results)
end
