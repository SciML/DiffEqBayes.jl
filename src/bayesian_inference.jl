export StanModel, bayesian_inference

struct StanModel{R,C}
  return_code::R
  chain_results::C
end
function bayesian_inference(prob::DEProblem,t,data;alg="integrate_ode_rk45",num_samples=1, num_warmup=1,kwargs...)
  length_of_y = string(length(prob.u0))
  length_of_parameter = string(length(prob.f.params))

  const parameter_estimation_model = "
  functions {
    real[] sho(real t,real[] u,real[] theta,real[] x_r,int[] x_i) {
      real du[$length_of_y];  // 2 = length(prob.u0) = length_of_y
      // placeholder for differential equation here
      du[1] = theta[1] * u[1] - 1.0*u[1] * u[2]; //string(prob.f.pfuncs[1])
      du[2] = -3.0*u[1] + 1.0*u[1] * u[2];
      return du;
      }
    }
  data {
    real u0[$length_of_y]; // 2 = length(prob.u0)
    int<lower=1> T;
    real u[T,$length_of_y]; // 2 = length(prob.u0)
    real t0;
    real ts[T];
  }
  transformed data {
    real x_r[0];
    int x_i[0];
  }
  parameters {
    vector<lower=0>[$length_of_y] sigma;   // 2 = length(prob.u0)
    real theta[$length_of_parameter];   // // 1=length(prob.f.params)
  }
  model{
    real u_hat[T,$length_of_y]; // 2 = length(prob.u0)
    sigma ~ inv_gamma(2, 3);
    // placeholder for priors here
    theta[1] ~ normal(1.5, 1); //1=length(prob.f.params), 1.5=prob.f.a
    u_hat = $alg(sho, u0, t0, ts, theta, x_r, x_i);
    for (t in 1:T){
      u[t] ~ normal(u_hat[t], sigma);
      }
  }
  "

  stanmodel = Stanmodel(num_samples=num_samples, num_warmup=num_warmup, name="parameter_estimation_model", model=parameter_estimation_model);
  const parameter_estimation_data = Dict("y0"=>prob.u0, "T" => size(t)[1], "y" => data', "t0" => prob.tspan[1], "ts"=>t)
  return_code, chain_results = stan(stanmodel, [parameter_estimation_data]; CmdStanDir=CMDSTAN_HOME)
  return StanModel(return_code,chain_results)
end
