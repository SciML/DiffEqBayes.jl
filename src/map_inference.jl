struct MAP
    prob
    alg
    priors
    t
    data_distributions
    kwargs
end

function (f::MAP)(p)
    prob = problem_new_parameters(f.prob, p)
    sol = solve(prob,f.alg;saveat=f.t)
    distributions = f.data_distributions
    fill_length = length(f.t)-length(sol)
    for i in 1:fill_length
      push!(sol.u,fill(Inf,size(sol[1])))
    end
    ll = 0.0
    if eltype(priors) <: UnivariateDistribution
        for j in 1:length(f.priors)
            ll -= logpdf(priors[j],f.prob.p[i])
        end
    else
        ll -= logpdf(priors,f.prob.p)
    end
      
    if eltype(distributions) <: UnivariateDistribution
        for j in 1:length(f.t), i in 1:length(sol[1][1])
            # i is the number of time points
            # j is the size of the system
            # corresponds to distributions[i,j]
            ll -= logpdf(distributions[i,j],sol[i,j])
        end
    else # MultivariateDistribution
        for j in 1:length(f.t), i in 1:length(sol[1][1])
            # i is the number of time points
            # j is the size of the system
            # corresponds to distributions[i,j]
            ll -= logpdf(distributions[i],sol[i])
        end
    end
      
    ll
end      

function map_inference(prob::ODEProblem,alg,t,data_distributions,priors; num_samples=1000, kwargs...)
    map_objective = MAP(prob,alg,priors,t,data_distributions,kwargs)
    map_objective
end

    
    
    